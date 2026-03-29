# NLP Tutorial Series

A series of 15 Jupyter notebooks taking you from "beginner in nonlinear programming" to understanding every component of ripopt — written in Python (numpy, scipy, matplotlib), with each notebook connecting directly to ripopt's source code.

All notebooks are in [`tutorials/`](https://github.com/jkitchin/ripopt/tree/main/tutorials) in the repository.

## Module 1: NLP Foundations

| # | Notebook | Key concepts |
|---|---|---|
| 01 | [Introduction to Nonlinear Programming](https://github.com/jkitchin/ripopt/blob/main/tutorials/01_nlp_intro.ipynb) | Problem formulation, feasibility vs optimality, HS071 |
| 02 | [Unconstrained Optimization](https://github.com/jkitchin/ripopt/blob/main/tutorials/02_unconstrained.ipynb) | Gradient descent, Newton's method, Armijo/Wolfe line search, L-BFGS |
| 03 | [KKT Conditions](https://github.com/jkitchin/ripopt/blob/main/tutorials/03_kkt_conditions.ipynb) | Lagrangian, complementarity, geometric interpretation, LICQ |

## Module 2: Classical Constrained Methods

| # | Notebook | Key concepts |
|---|---|---|
| 04 | [Penalty Methods and Augmented Lagrangian](https://github.com/jkitchin/ripopt/blob/main/tutorials/04_penalty_augmented_lagrangian.ipynb) | Quadratic penalty, ill-conditioning, Method of Multipliers |
| 05 | [Sequential Quadratic Programming](https://github.com/jkitchin/ripopt/blob/main/tutorials/05_sqp.ipynb) | QP subproblem, active set, BFGS update, L1 merit function |

## Module 3: Barrier Methods

| # | Notebook | Key concepts |
|---|---|---|
| 06 | [The Logarithmic Barrier Method](https://github.com/jkitchin/ripopt/blob/main/tutorials/06_barrier_methods.ipynb) | Log barrier, central path, barrier parameter sequence |
| 07 | [Primal-Dual Interior Point Method](https://github.com/jkitchin/ripopt/blob/main/tutorials/07_primal_dual_ipm.ipynb) | Perturbed KKT, Newton step derivation, fraction-to-boundary, complete IPM |

## Module 4: Linear Algebra at the Core

| # | Notebook | Key concepts |
|---|---|---|
| 08 | [The KKT (Saddle-Point) Matrix](https://github.com/jkitchin/ripopt/blob/main/tutorials/08_kkt_matrix.ipynb) | KKT structure, Sylvester's Law of Inertia, condensed KKT |
| 09 | [LDL^T and Bunch-Kaufman Pivoting](https://github.com/jkitchin/ripopt/blob/main/tutorials/09_ldlt_bunch_kaufman.ipynb) | Symmetric Gaussian elimination, BK pivoting, inertia from D |
| 10 | [Inertia Correction and Regularization](https://github.com/jkitchin/ripopt/blob/main/tutorials/10_inertia_correction.ipynb) | δ_w and δ_c perturbations, correction loop, iterative refinement |

## Module 5: Line Search, Filter, and Convergence

| # | Notebook | Key concepts |
|---|---|---|
| 11 | [Filter Line Search](https://github.com/jkitchin/ripopt/blob/main/tutorials/11_filter_line_search.ipynb) | Filter concept, acceptance criteria, switching condition, SOC |
| 12 | [Convergence Criteria and Scaling](https://github.com/jkitchin/ripopt/blob/main/tutorials/12_convergence_criteria.ipynb) | Three KKT residuals, s_d/s_c scaling, NLP scaling |

## Module 6: Robustness

| # | Notebook | Key concepts |
|---|---|---|
| 13 | [The Restoration Phase](https://github.com/jkitchin/ripopt/blob/main/tutorials/13_restoration.ipynb) | Gauss-Newton restoration, NLP restoration subproblem, recovery strategies |
| 14 | [Mehrotra Predictor-Corrector](https://github.com/jkitchin/ripopt/blob/main/tutorials/14_mehrotra_gondzio.ipynb) | Affine-scaling predictor, adaptive centering σ = (μ_aff/μ)³, Gondzio corrections |

## Module 7: ripopt in Practice

| # | Notebook | Key concepts |
|---|---|---|
| 15 | [ripopt in Practice](https://github.com/jkitchin/ripopt/blob/main/tutorials/15_ripopt_in_practice.ipynb) | NlpProblem interface, fallback cascade, sensitivity analysis, benchmarking |

## Running the notebooks

```bash
git clone https://github.com/jkitchin/ripopt.git
cd ripopt/tutorials
pip install numpy scipy matplotlib jupyter
jupyter lab
```

Each notebook is self-contained. No Rust installation required for notebooks 01–14. Notebook 15 describes the Rust interface with Python equivalents for all concepts.

## ripopt source pointers

Each notebook ends with a "Connection to ripopt" section. The key mappings:

| Topic | Notebook | ripopt source |
|---|---|---|
| KKT assembly | 08 | `src/kkt.rs` — `assemble_kkt()` |
| Bunch-Kaufman LDL^T | 09 | `src/linear_solver/dense.rs` — `DenseLdl` |
| Inertia correction | 10 | `src/kkt.rs` — `factor_with_inertia_correction()` |
| Filter line search | 11 | `src/filter.rs` — `check_acceptability()` |
| Convergence | 12 | `src/convergence.rs` — `check_convergence()` |
| Restoration | 13 | `src/restoration.rs`, `src/restoration_nlp.rs` |
| Mehrotra PC | 14 | `src/ipm.rs` — main loop |
| Problem interface | 15 | `src/problem.rs` — `NlpProblem` trait |
