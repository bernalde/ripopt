/**
 * gams_ripopt.c — GAMS solver link for ripopt.
 *
 * Translates between the GAMS Modeling Object (GMO) API and ripopt's C API
 * (which mirrors Ipopt). Produces a shared library that drops into a GAMS
 * installation so that models can use `option nlp = ripopt;`.
 *
 * Entry points (prefix "rip" registered in gmscmpun.txt):
 *   ripCreate      — allocate solver data
 *   ripFree        — free solver data
 *   ripReadyAPI    — initialize GAMS API libraries
 *   ripCallSolver  — extract problem, solve, report solution
 *
 * Build:  make -C gams  (see Makefile)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#ifndef _WIN32
#include <dlfcn.h>
#endif

#include "ripopt.h"
#include "gmomcc.h"
#include "gevmcc.h"

/* ---------------------------------------------------------------------------
 * Solver data
 * --------------------------------------------------------------------------- */

typedef struct {
    gmoHandle_t gmo;
    gevHandle_t gev;

    /* Problem dimensions */
    int n;          /* number of variables */
    int m;          /* number of constraints */
    int nnz_jac;    /* Jacobian nonzeros */
    int nnz_hess;   /* Hessian nonzeros (lower triangle) */

    /* 1.0 for minimization, -1.0 for maximization */
    double obj_sign;

    /* Jacobian CSR structure (from gmoGetMatrixRow) */
    int *jac_rowstart;  /* length m + 1 */
    int *jac_colidx;    /* length nnz_jac */

    /* Dense gradient scratch buffer (length n) */
    double *grad_buf;

    /* Has analytical Hessian? */
    int have_hessian;
} RipoptGamsData;

/* ---------------------------------------------------------------------------
 * Callback: objective f(x)
 * --------------------------------------------------------------------------- */

static int gams_eval_f(int n, const double *x, int new_x,
                       double *obj, void *user_data)
{
    RipoptGamsData *d = (RipoptGamsData *)user_data;
    double fval;
    int numerr;

    (void)n; (void)new_x;

    if (gmoEvalFuncObj(d->gmo, x, &fval, &numerr) != 0 || numerr > 0)
        return 0;

    *obj = d->obj_sign * fval;
    return 1;
}

/* ---------------------------------------------------------------------------
 * Callback: gradient of f(x)
 * --------------------------------------------------------------------------- */

static int gams_eval_grad_f(int n, const double *x, int new_x,
                            double *grad, void *user_data)
{
    RipoptGamsData *d = (RipoptGamsData *)user_data;
    double fval, gx;
    int numerr;

    (void)new_x;

    if (gmoEvalGradObj(d->gmo, x, &fval, grad, &gx, &numerr) != 0 || numerr > 0)
        return 0;

    if (d->obj_sign < 0.0) {
        for (int j = 0; j < n; j++)
            grad[j] = -grad[j];
    }
    return 1;
}

/* ---------------------------------------------------------------------------
 * Callback: constraints g(x)
 * --------------------------------------------------------------------------- */

static int gams_eval_g(int n, const double *x, int new_x,
                       int m, double *g, void *user_data)
{
    RipoptGamsData *d = (RipoptGamsData *)user_data;
    double fval;
    int numerr;

    (void)n; (void)new_x;

    for (int i = 0; i < m; i++) {
        if (gmoEvalFunc(d->gmo, i, x, &fval, &numerr) != 0 || numerr > 0)
            return 0;
        g[i] = fval;
    }
    return 1;
}

/* ---------------------------------------------------------------------------
 * Callback: Jacobian of constraints
 *
 * Structure mode (values == NULL): expand CSR to COO
 * Values mode: evaluate each row's dense gradient, extract sparse entries
 * --------------------------------------------------------------------------- */

static int gams_eval_jac_g(int n, const double *x, int new_x,
                           int m, int nele_jac,
                           int *iRow, int *jCol,
                           double *values, void *user_data)
{
    RipoptGamsData *d = (RipoptGamsData *)user_data;

    (void)n; (void)new_x; (void)nele_jac;

    if (values == NULL) {
        /* Sparsity pattern: expand CSR to COO (0-based) */
        int k = 0;
        for (int i = 0; i < m; i++) {
            for (int j = d->jac_rowstart[i]; j < d->jac_rowstart[i + 1]; j++) {
                iRow[k] = i;
                jCol[k] = d->jac_colidx[j];
                k++;
            }
        }
    } else {
        /* Values: evaluate gradient row by row */
        int k = 0;
        for (int i = 0; i < m; i++) {
            double fval, gx;
            int numerr;
            memset(d->grad_buf, 0, d->n * sizeof(double));

            if (gmoEvalGrad(d->gmo, i, x, &fval, d->grad_buf, &gx, &numerr) != 0
                || numerr > 0)
                return 0;

            for (int j = d->jac_rowstart[i]; j < d->jac_rowstart[i + 1]; j++) {
                values[k] = d->grad_buf[d->jac_colidx[j]];
                k++;
            }
        }
    }
    return 1;
}

/* ---------------------------------------------------------------------------
 * Callback: Hessian of the Lagrangian (lower triangle)
 *
 * ripopt convention:  H = obj_factor * nabla^2 f + sum_i lambda_i * nabla^2 c_i
 *
 * gmoHessLagValue(gmo, x, pi, w, objweight, conweight, numerr) computes:
 *   w = objweight * nabla^2 f + conweight * sum_i pi_i * nabla^2 c_i
 *
 * GAMS multiplier convention: pi_gams = -lambda_ipopt.
 * The Ipopt solver link negates lambda before passing to gmoHessLagValue.
 * Equivalently, we pass lambda directly with conweight = -1.0:
 *   w = objweight * nabla^2 f + (-1) * sum_i lambda_i * nabla^2 c_i
 *     = objweight * nabla^2 f - sum_i lambda_i * nabla^2 c_i
 * which matches negating lambda and using conweight = 1.0.
 *
 * For maximization, obj_sign = -1, so objweight = -obj_factor negates the
 * objective Hessian (ripopt minimizes -f).
 * --------------------------------------------------------------------------- */

static int gams_eval_h(int n, const double *x, int new_x,
                       double obj_factor,
                       int m, const double *lambda, int new_lambda,
                       int nele_hess,
                       int *iRow, int *jCol,
                       double *values, void *user_data)
{
    RipoptGamsData *d = (RipoptGamsData *)user_data;

    (void)n; (void)new_x; (void)m; (void)new_lambda; (void)nele_hess;

    if (values == NULL) {
        /* Sparsity pattern: COO (lower triangle) from gmoHessLagStruct */
        gmoHessLagStruct(d->gmo, iRow, jCol);
    } else {
        int numerr;
        double objweight = d->obj_sign * obj_factor;

        /* conweight = -1.0: equivalent to negating lambda (GAMS sign convention) */
        if (gmoHessLagValue(d->gmo, x, lambda, values,
                            objweight, -1.0, &numerr) != 0
            || numerr > 0)
            return 0;
    }
    return 1;
}

/* ---------------------------------------------------------------------------
 * Option file parsing
 *
 * Reads lines of the form "key value" from ripopt.opt (or .op2 etc.)
 * Lines starting with '*' are comments. Blank lines are skipped.
 * --------------------------------------------------------------------------- */

static void parse_option_file(RipoptProblem nlp, const char *filename,
                              gevHandle_t gev)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) return;

    char line[512];
    while (fgets(line, sizeof(line), fp)) {
        /* Skip comments and blank lines */
        char *p = line;
        while (*p && isspace((unsigned char)*p)) p++;
        if (*p == '\0' || *p == '*' || *p == '#') continue;

        /* Remove trailing newline */
        char *nl = strchr(p, '\n');
        if (nl) *nl = '\0';
        char *cr = strchr(p, '\r');
        if (cr) *cr = '\0';

        /* Split "key value" */
        char key[256], val[256];
        if (sscanf(p, "%255s %255s", key, val) < 2)
            continue;

        /* Try integer first, then double, then string */
        char *endptr;
        long ival = strtol(val, &endptr, 10);
        if (*endptr == '\0') {
            if (ripopt_add_int_option(nlp, key, (int)ival)) {
                gevLogStat(gev, line);
                continue;
            }
        }

        double dval = strtod(val, &endptr);
        if (*endptr == '\0') {
            if (ripopt_add_num_option(nlp, key, dval)) {
                gevLogStat(gev, line);
                continue;
            }
        }

        if (ripopt_add_str_option(nlp, key, val)) {
            gevLogStat(gev, line);
            continue;
        }

        /* Unknown option */
        char msgbuf[512];
        snprintf(msgbuf, sizeof(msgbuf), "*** Warning: unknown option '%s'", key);
        gevLogStat(gev, msgbuf);
    }
    fclose(fp);
}

/* ---------------------------------------------------------------------------
 * GAMS solver link entry points
 * --------------------------------------------------------------------------- */

#if defined(_WIN32)
# define DllExport __declspec(dllexport)
# define STDCALL   __stdcall
#else
# define DllExport __attribute__((visibility("default")))
# define STDCALL
#endif

/** Allocate solver data and initialize GAMS API wrappers.
 *
 * Returns 0 on success, 1 on failure (with error in msgBuf).
 * GAMS calls gmoGetReady/gevGetReady here (before ReadyAPI),
 * matching the pattern used by IPOPT and other solver links.
 */
DllExport int STDCALL ripCreate(void **Cptr, char *msgBuf, int msgBufLen)
{
    *Cptr = NULL;

    /* Initialize GAMS API wrappers (function pointers) */
    if (!gmoGetReady(msgBuf, msgBufLen))
        return 1;
    if (!gevGetReady(msgBuf, msgBufLen))
        return 1;

    RipoptGamsData *data = (RipoptGamsData *)calloc(1, sizeof(RipoptGamsData));
    if (!data) {
        snprintf(msgBuf, msgBufLen, "ripopt: memory allocation failed");
        return 1;
    }
    *Cptr = data;
    msgBuf[0] = '\0';
    return 0;
}

/** Free solver data. */
DllExport void STDCALL ripFree(void **Cptr)
{
    if (Cptr && *Cptr) {
        RipoptGamsData *data = (RipoptGamsData *)*Cptr;
        free(data->jac_rowstart);
        free(data->jac_colidx);
        free(data->grad_buf);
        free(data);
        *Cptr = NULL;
    }
}

/** Initialize GAMS API — receive GMO handle and extract GEV. */
DllExport int STDCALL ripReadyAPI(void *Cptr, gmoHandle_t gmo)
{
    RipoptGamsData *data = (RipoptGamsData *)Cptr;

    data->gmo = gmo;
    data->gev = (gevHandle_t)gmoEnvironment(gmo);
    if (!data->gev) {
        gmoSolveStatSet(gmo, gmoSolveStat_SetupErr);
        gmoModelStatSet(gmo, gmoModelStat_ErrorNoSolution);
        return 1;
    }

    return 0;
}

/** Map ripopt return status to GAMS model/solve status. */
static void map_status_to_gams(int ripopt_status, int *model_stat, int *solve_stat)
{
    switch (ripopt_status) {
    case 0:  /* SOLVE_SUCCEEDED */
        *model_stat = gmoModelStat_OptimalLocal;
        *solve_stat = gmoSolveStat_Normal;
        break;
    case 1:  /* ACCEPTABLE_LEVEL */
        *model_stat = gmoModelStat_Feasible;
        *solve_stat = gmoSolveStat_Normal;
        break;
    case 2:  /* INFEASIBLE_PROBLEM */
        *model_stat = gmoModelStat_InfeasibleLocal;
        *solve_stat = gmoSolveStat_Normal;
        break;
    case 5:  /* MAXITER_EXCEEDED */
        *model_stat = gmoModelStat_InfeasibleIntermed;
        *solve_stat = gmoSolveStat_Iteration;
        break;
    case 6:  /* RESTORATION_FAILED */
        *model_stat = gmoModelStat_InfeasibleLocal;
        *solve_stat = gmoSolveStat_Solver;
        break;
    case 7:  /* ERROR_IN_STEP */
        *model_stat = gmoModelStat_ErrorNoSolution;
        *solve_stat = gmoSolveStat_EvalError;
        break;
    case 10: /* NOT_ENOUGH_DOF */
    case 11: /* INVALID_PROBLEM */
        *model_stat = gmoModelStat_ErrorNoSolution;
        *solve_stat = gmoSolveStat_SetupErr;
        break;
    default: /* INTERNAL_ERROR or unknown */
        *model_stat = gmoModelStat_ErrorNoSolution;
        *solve_stat = gmoSolveStat_InternalErr;
        break;
    }
}

/** Solve the NLP problem. */
DllExport int STDCALL ripCallSolver(void *Cptr)
{
    RipoptGamsData *data = (RipoptGamsData *)Cptr;
    gmoHandle_t gmo = data->gmo;
    gevHandle_t gev = data->gev;
    char msg[512];
    int n, m, rc;

    /* All heap pointers initialized here so that goto cleanup is always safe */
    double *x_l     = NULL;
    double *x_u     = NULL;
    double *g_l     = NULL;
    double *g_u     = NULL;
    double *x       = NULL;
    double *g_vals  = NULL;
    double *mult_g  = NULL;
    double *mult_xl = NULL;
    double *mult_xu = NULL;
    RipoptProblem nlp = NULL;

    gevLogStat(gev, "");
    gevLogStat(gev, "--- RIPOPT: A Rust Interior-Point Optimizer");
    gevLogStat(gev, "");

    /* ---------------------------------------------------------------
     * Configure GMO
     * --------------------------------------------------------------- */
    gmoObjStyleSet(gmo, gmoObjType_Fun);
    gmoObjReformSet(gmo, 1);
    gmoIndexBaseSet(gmo, 0);

    /* Objective sense */
    data->obj_sign = (gmoSense(gmo) == gmoObj_Max) ? -1.0 : 1.0;

    /* ---------------------------------------------------------------
     * Extract problem dimensions
     * --------------------------------------------------------------- */
    data->n = gmoN(gmo);
    data->m = gmoM(gmo);
    data->nnz_jac = gmoNZ(gmo);
    n = data->n;
    m = data->m;

    if (n == 0) {
        gevLogStat(gev, "*** Error: problem has no variables");
        gmoSolveStatSet(gmo, gmoSolveStat_SetupErr);
        gmoModelStatSet(gmo, gmoModelStat_ErrorNoSolution);
        return 1;
    }

    snprintf(msg, sizeof(msg), "  Variables: %d, Constraints: %d, Jacobian NZ: %d",
             n, m, data->nnz_jac);
    gevLogStat(gev, msg);

    /* ---------------------------------------------------------------
     * Load Hessian
     * --------------------------------------------------------------- */
    {
        int do2dir, doHess;
        rc = gmoHessLoad(gmo, 0.0, &do2dir, &doHess);
        if (rc != 0 || !doHess) {
            data->have_hessian = 0;
            data->nnz_hess = 0;
            gevLogStat(gev, "  Analytical Hessian not available, using L-BFGS");
        } else {
            data->have_hessian = 1;
            data->nnz_hess = gmoHessLagNz(gmo);
            snprintf(msg, sizeof(msg), "  Hessian NZ: %d", data->nnz_hess);
            gevLogStat(gev, msg);
        }
    }

    /* ---------------------------------------------------------------
     * Variable bounds
     * --------------------------------------------------------------- */
    x_l = (double *)malloc(n * sizeof(double));
    x_u = (double *)malloc(n * sizeof(double));
    if (!x_l || !x_u) goto oom;

    gmoGetVarLower(gmo, x_l);
    gmoGetVarUpper(gmo, x_u);

    /* Map GAMS infinity to ripopt infinity (1e30) */
    {
        double gams_pinf = gmoPinf(gmo);
        double gams_minf = gmoMinf(gmo);
        for (int j = 0; j < n; j++) {
            if (x_l[j] <= gams_minf) x_l[j] = -1e30;
            if (x_u[j] >= gams_pinf) x_u[j] =  1e30;
        }
    }

    /* ---------------------------------------------------------------
     * Constraint bounds from equation types and RHS
     * --------------------------------------------------------------- */
    if (m > 0) {
        g_l = (double *)malloc(m * sizeof(double));
        g_u = (double *)malloc(m * sizeof(double));
        if (!g_l || !g_u) goto oom;

        for (int i = 0; i < m; i++) {
            int etyp = gmoGetEquTypeOne(gmo, i);
            double rhs = gmoGetRhsOne(gmo, i);

            switch (etyp) {
            case gmoequ_E:  /* =E= equality */
                g_l[i] = rhs;
                g_u[i] = rhs;
                break;
            case gmoequ_G:  /* =G= greater-or-equal */
                g_l[i] = rhs;
                g_u[i] = 1e30;
                break;
            case gmoequ_L:  /* =L= less-or-equal */
                g_l[i] = -1e30;
                g_u[i] = rhs;
                break;
            case gmoequ_N:  /* =N= free / nonbinding */
                g_l[i] = -1e30;
                g_u[i] = 1e30;
                break;
            default:
                snprintf(msg, sizeof(msg),
                         "*** Warning: unsupported equation type %d for row %d",
                         etyp, i);
                gevLogStat(gev, msg);
                g_l[i] = -1e30;
                g_u[i] = 1e30;
                break;
            }
        }
    }

    /* ---------------------------------------------------------------
     * Jacobian structure (CSR from GMO, stored for value callbacks)
     * --------------------------------------------------------------- */
    if (m > 0 && data->nnz_jac > 0) {
        data->jac_rowstart = (int *)malloc((m + 1) * sizeof(int));
        data->jac_colidx   = (int *)malloc(data->nnz_jac * sizeof(int));
        data->grad_buf     = (double *)calloc(n, sizeof(double));
        if (!data->jac_rowstart || !data->jac_colidx || !data->grad_buf)
            goto oom;

        /* Get CSR structure; jacval and nlflag are unused after setup */
        double *jacval_tmp = (double *)malloc(data->nnz_jac * sizeof(double));
        int    *nlflag_tmp = (int *)malloc(data->nnz_jac * sizeof(int));
        if (!jacval_tmp || !nlflag_tmp) {
            free(jacval_tmp);
            free(nlflag_tmp);
            goto oom;
        }

        gmoGetMatrixRow(gmo, data->jac_rowstart, data->jac_colidx,
                        jacval_tmp, nlflag_tmp);
        free(jacval_tmp);
        free(nlflag_tmp);
    }

    /* ---------------------------------------------------------------
     * Create ripopt problem
     * --------------------------------------------------------------- */
    {
        Eval_H_CB hess_cb = data->have_hessian ? gams_eval_h : NULL;
        int nnz_hess_arg  = data->have_hessian ? data->nnz_hess : 0;

        nlp = ripopt_create(
            n, x_l, x_u,
            m, g_l, g_u,
            data->nnz_jac, nnz_hess_arg,
            gams_eval_f,
            gams_eval_grad_f,
            gams_eval_g,
            gams_eval_jac_g,
            hess_cb);
    }

    if (!nlp) {
        gevLogStat(gev, "*** Error: ripopt_create failed");
        gmoSolveStatSet(gmo, gmoSolveStat_SetupErr);
        gmoModelStatSet(gmo, gmoModelStat_ErrorNoSolution);
        rc = 1;
        goto cleanup;
    }

    /* ---------------------------------------------------------------
     * Default options from GAMS environment
     * --------------------------------------------------------------- */
    {
        int iterlim = gevGetIntOpt(gev, gevIterLim);
        if (iterlim < ITERLIM_INFINITY)
            ripopt_add_int_option(nlp, "max_iter", iterlim);

        double reslim = gevGetDblOpt(gev, gevResLim);
        if (reslim < RESLIM_INFINITY)
            ripopt_add_num_option(nlp, "max_wall_time", reslim);
    }

    /* Default print level */
    ripopt_add_int_option(nlp, "print_level", 5);

    /* Use L-BFGS if no analytical Hessian */
    if (!data->have_hessian)
        ripopt_add_str_option(nlp, "hessian_approximation", "limited-memory");

    /* ---------------------------------------------------------------
     * Read option file (ripopt.opt, ripopt.op2, ...)
     * --------------------------------------------------------------- */
    if (gmoOptFile(gmo) > 0) {
        char optfilename[512];
        gmoNameOptFile(gmo, optfilename);
        snprintf(msg, sizeof(msg), "  Reading option file %s", optfilename);
        gevLogStat(gev, msg);
        parse_option_file(nlp, optfilename, gev);
    }

    /* ---------------------------------------------------------------
     * Allocate solution arrays and set initial point
     * --------------------------------------------------------------- */
    x       = (double *)malloc(n * sizeof(double));
    g_vals  = m > 0 ? (double *)malloc(m * sizeof(double)) : NULL;
    mult_g  = m > 0 ? (double *)calloc(m, sizeof(double)) : NULL;
    mult_xl = (double *)calloc(n, sizeof(double));
    mult_xu = (double *)calloc(n, sizeof(double));

    if (!x || (m > 0 && (!g_vals || !mult_g)) || !mult_xl || !mult_xu)
        goto oom;

    gmoGetVarL(gmo, x);

    /* ---------------------------------------------------------------
     * Solve
     * --------------------------------------------------------------- */
    {
        double obj_val = 0.0;
        int status = ripopt_solve(nlp, x, g_vals, &obj_val,
                                  mult_g, mult_xl, mult_xu,
                                  (void *)data);

        /* -----------------------------------------------------------
         * Report solution back to GAMS
         * ----------------------------------------------------------- */
        int model_stat, solve_stat;
        map_status_to_gams(status, &model_stat, &solve_stat);

        /* Objective in GAMS convention (undo our sign flip for max) */
        if (status == 0 || status == 1) {
            double gams_obj = (data->obj_sign < 0.0) ? -obj_val : obj_val;
            gmoSetHeadnTail(gmo, gmoHobjval, gams_obj);
        }

        gmoModelStatSet(gmo, model_stat);
        gmoSolveStatSet(gmo, solve_stat);

        /* Negate constraint multipliers: ripopt lambda -> GAMS pi */
        if (mult_g) {
            for (int i = 0; i < m; i++)
                mult_g[i] = -mult_g[i];
        }

        /* Set primal + dual solution */
        gmoSetSolution2(gmo, x, mult_g);

        /* Variable marginals: z_L - z_U, negated for max problems */
        {
            double *var_marg = (double *)calloc(n, sizeof(double));
            if (var_marg) {
                for (int j = 0; j < n; j++)
                    var_marg[j] = mult_xl[j] - mult_xu[j];
                if (data->obj_sign < 0.0) {
                    for (int j = 0; j < n; j++)
                        var_marg[j] = -var_marg[j];
                }
                gmoSetVarM(gmo, var_marg);
                free(var_marg);
            }
        }

        snprintf(msg, sizeof(msg),
                 "  Solve status: %d (%s), Model status: %d (%s)",
                 solve_stat,
                 solve_stat <= 13 ? solveStatusTxt[solve_stat] : "?",
                 model_stat,
                 model_stat <= 19 ? modelStatusTxt[model_stat] : "?");
        gevLogStat(gev, msg);
    }

    rc = 0;

cleanup:
    if (nlp) ripopt_free(nlp);
    free(x_l);
    free(x_u);
    free(g_l);
    free(g_u);
    free(x);
    free(g_vals);
    free(mult_g);
    free(mult_xl);
    free(mult_xu);
    if (data->have_hessian)
        gmoHessUnload(gmo);
    return rc;

oom:
    gevLogStat(gev, "*** Error: out of memory");
    gmoSolveStatSet(gmo, gmoSolveStat_InternalErr);
    gmoModelStatSet(gmo, gmoModelStat_ErrorNoSolution);
    rc = 1;
    goto cleanup;
}
