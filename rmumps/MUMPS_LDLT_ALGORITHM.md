Good. Now I have everything I need. Let me compose the complete specification.

---

Here is the complete algorithmic specification of MUMPS's frontal LDLT factorization for a single type-1 supernode, traced directly from the source code.

## Notation and Setup

All source references are to files under `/Users/jkitchin/Dropbox/projects/ripopt/ref/mumps/src/`. The three core subroutines are:

- `DMUMPS_FAC_I_LDLT` in `dfac_front_aux.F` lines 1147-1676 — pivot search and swap
- `DMUMPS_FAC_MQ_LDLT` in `dfac_front_aux.F` lines 1677-1989 — within-panel rank-1/rank-2 update
- `DMUMPS_FAC_SQ_LDLT` in `dfac_front_aux.F` lines 1990-2089 — between-panel BLAS-3 update
- `DMUMPS_SWAP_LDLT` in `dfac_front_aux.F` lines 2090-2151 — symmetric row/column swap
- `DMUMPS_FAC_LDLT_COPY2U_SCALEL` in `dfac_front_aux.F` lines 2152-2265 — copy L to U and scale by D^{-1}

The driver loop is in `DMUMPS_FAC1_LDLT` in `dfac_front_LDLT_type1.F`.

### Index convention

MUMPS uses **1-based** Fortran indexing throughout. The frontal matrix `A` is stored as a flat 1D array starting at position `POSELT`. For a matrix element at row `r`, column `c` (both 1-based), the position in `A` is:

```
A[POSELT + (c-1)*LDA + (r-1)]
```

where `LDA = NFRONT`. The matrix is **column-major**, full dense square of size `NFRONT x NFRONT`.

### Key variables

- `NPIV` — number of pivots already eliminated (read from `IW(IOLDPS+1+XSIZE)`)
- `NPIVP1 = NPIV + 1` — the position where the next pivot will go
- `NASS` — number of fully-summed (assemblable) columns
- `NFRONT` — total frontal matrix size
- `NCB = NFRONT - NASS` — contribution block size (rows/cols not eliminated at this node)
- `POSELT` — offset of (1,1) element in A array

### Frontal matrix regions

```
         col 1..NPIV     NPIV+1..NASS    NASS+1..NFRONT
row 1..NPIV    [already factored]
NPIV+1..NASS   [L entries | U saved]   [active FS block]   [FS-CB coupling]
NASS+1..NFRONT [L entries | U saved]   [CB-FS coupling]    [contribution block]
```

Columns 1..NPIV are already factored. The "active" region starts at row NPIV+1, column NPIV+1.

---

## PHASE 0: LOOP STRUCTURE (The Driver)

Source: `dfac_front_LDLT_type1.F` lines 150-783.

For non-BLR, the structure simplifies to:

```
NBKJIB_ORIG = inner_block_size(NASS)  // typically 16 or 32
NBLR_ORIG = KEEP(420)                 // typically 128 (= 4*KEEP(6))

NPIV = 0
IEND_BLR = 0          // end of current outer panel
IEND_BLOCK = 0        // end of current inner block

// === OUTER PANEL LOOP ===
while IEND_BLR < NASS:
    IEND_BLR = min(IEND_BLR + NBLR_ORIG, NASS)

    // === INNER BLOCK LOOP ===
    while IEND_BLOCK < IEND_BLR:
        IBEG_BLOCK = NPIV + 1     // read from IW(IOLDPS+1+XSIZE) + 1
        IEND_BLOCK = min(IEND_BLOCK + NBKJIB_ORIG, IEND_BLR)

        // === PIVOT-BY-PIVOT LOOP within this inner block ===
        loop:
            // STEP A: Find pivot
            call FAC_I_LDLT(...)  -> INOPV, PIVSIZ

            if INOPV == 1:  // no pivot found, search exhausted NASS
                if STATICMODE:
                    INOPV = -1
                    goto STEP_A   // retry with static pivoting
                else:
                    LASTPANEL = true
                    break inner block loop AND outer panel loop
            else if INOPV == 2:  // no pivot in this inner block, but more blocks
                break pivot loop (to advance to next inner block)
            // INOPV <= 0: pivot found (0 = normal, -1 = static)

            // STEP B: Rank-1 or rank-2 update within this inner block
            LAST_ROW = NFRONT  // for PIVOT_OPTION=3
            NVSCHUR_K253 = NVSCHUR + KEEP(253)
            call FAC_MQ_LDLT(...)  -> IFINB

            // Mark 2x2 pivot in IW: negate column index of second pivot column
            if PIVSIZ == 2:
                IW(IOLDPS + NPIV + 7 + XSIZE + NFRONT) = -IW(...)

            NPIV = NPIV + PIVSIZ   // stored in IW(IOLDPS+1+XSIZE)

            if IFINB == 0:
                continue pivot loop (more columns in this inner block)
            else if IFINB == 1:
                break pivot loop (inner block exhausted, but panel continues)
            else if IFINB == -1:
                LASTPANEL = true
                break inner block loop (panel exhausted = NASS reached)

        // STEP C: Between-panel BLAS-3 update
        if IEND_BLR > IEND_BLOCK:
            // Update remaining columns in this BLR panel
            call FAC_SQ_LDLT(IBEG_BLOCK, IEND_BLOCK, NPIV,
                             ..., IEND_BLR, LAST_ROW, ...)
    end inner block loop

    // STEP D: Inter-panel BLAS-3 update
    // Update columns from IEND_BLR+1 to NASS, and rows from NASS+1 to NFRONT
    call FAC_SQ_LDLT(IBEG_BLR, IEND_BLR, NPIV,
                     ..., NASS, NASS, NASS, LAST_ROW, ...)
end outer panel loop
```

The key insight: MUMPS uses a **two-level blocking** even without BLR. The outer panel (NBLR_ORIG ~ 128) determines the inter-panel BLAS-3 granularity. The inner block (NBKJIB_ORIG ~ 16-32) is the granularity for rank-1/rank-2 updates within a panel.

---

## PHASE 1: PIVOT SEARCH — FAC_I_LDLT

Source: `dfac_front_aux.F` lines 1147-1676.

### Inputs

- `NPIV` — current pivot count (read from `IW(IOLDPS+1+XSIZE)`)
- `IBEG_BLOCK, IEND_BLOCK` — range of the current inner block (1-based column indices)
- `PIVOT_OPTION = 3` — search rows/columns across entire front
- `UU = CNTL(1) = 0.01` — threshold
- `SEUIL = CNTL(4)` or computed from DKEEP(1) — tiny pivot threshold

### Static pivot mode (INOPV == -1 on entry)

If entering with `INOPV == -1`, this means the previous normal search failed and we are retrying with static pivoting. The code takes the diagonal element at position `(NPIV+1, NPIV+1)` unconditionally:

```
APOS = POSELT + NPIV*(LDA+1)   // diagonal of column NPIV+1
if |A[APOS]| < SEUIL:
    if A[APOS] >= 0:  A[APOS] = CSEUIL
    else:              A[APOS] = -CSEUIL; NNEG++
    // (CSEUIL = max(DKEEP(1), SEUIL))
else:
    // accept as-is
goto end (pivot accepted, PIVSIZ = 1)
```

### Normal pivot search (INOPV == 0 on entry)

```
PIVSIZ = 1
NPIVP1 = NPIV + 1

// Inextpiv acceleration: start scanning from last-known good position
ISHIFT = 0
if KEEP(206) >= 1:
    if Inextpiv > NPIVP1 and Inextpiv <= IEND_BLOCK:
        ISHIFT = Inextpiv - NPIVP1
    IPIV_END = IEND_BLOCK + ISHIFT   // wraparound search range
    // Quick check: if the diagonal at NPIVP1 already passes, skip shift
    if ISHIFT > 0 and IS_MAXFROMM_AVAIL:
        check if diagonal(NPIVP1) >= UU * MAXFROMM and > SEUIL
        if yes: ISHIFT = 0   // no need to skip ahead
    if ISHIFT > 0:
        IS_MAXFROMM_AVAIL = false

// === MAIN CANDIDATE LOOP ===
for IPIV_SHIFT = NPIVP1+ISHIFT to IPIV_END:
    // Wraparound: scan [Inextpiv..IEND_BLOCK] then [NPIVP1..Inextpiv-1]
    if IPIV_SHIFT <= IEND_BLOCK:
        IPIV = IPIV_SHIFT
    else:
        IPIV = IPIV_SHIFT - IEND_BLOCK - 1 + NPIVP1
        if IBEG_BLOCK == NPIVP1: break  // full wrap done

    // Position in matrix: row NPIVP1..LDA of column IPIV
    // More precisely: row of candidate = IPIV, column being examined = NPIV+1
    // APOS = start of column IPIV in the active row range
    APOS = POSELT + (IPIV-1)*LDA + NPIV     // A[NPIV+1, IPIV] 0-based
    POSPV1 = APOS + (IPIV - NPIVP1)         // A[IPIV, IPIV] = diagonal
    PIVOT = A[POSPV1]

    // ----- THRESHOLD=0 PATH (no pivoting) -----
    if UU == 0:
        // Accept diagonal unconditionally (with SEUIL replacement if tiny)
        if |PIVOT| < SEUIL: replace with ±CSEUIL
        goto accept

    // ----- MAXFROMM FAST-ACCEPT (from previous MQ update) -----
    if IS_MAXFROMM_AVAIL:
        // MAXFROMM is the max in the first trailing column from the MQ update
        if MAXFROMM > PIVNUL:
            if |PIVOT| >= UU * MAXFROMM and |PIVOT| > SEUIL:
                goto accept_1x1

    // ----- COMPUTE COLUMN MAXIMUM (AMAX, JMAX) -----
    // AMAX = max |A[k, IPIV]| for k in NPIVP1..IEND_BLOCK, k != IPIV
    // (i.e., the active fully-summed part of column IPIV)
    // By symmetry, column IPIV row range = A[NPIVP1..NFRONT, IPIV]
    // But below diagonal, A[IPIV, k] for k > IPIV is stored in column k
    
    // Scan within the active block [NPIVP1..IEND_BLOCK]:
    AMAX = -1
    JMAX = 0

    // Part 1: rows NPIVP1 to IPIV-1 (below diagonal of column IPIV,
    //          stored in row IPIV of columns NPIVP1..IPIV-1)
    for JJ = APOS to POSPV1-1:    // these are A[NPIVP1..IPIV-1, IPIV]
        if |A[JJ]| > AMAX:
            AMAX = |A[JJ]|
            JMAX = IPIV - (POSPV1 - JJ)  // column index of max

    // Part 2: rows IPIV+1 to IEND_BLOCK (stored in columns IPIV+1..IEND_BLOCK,
    //          at row IPIV)
    J1 = POSPV1 + LDA
    for J = 1 to IEND_BLOCK - IPIV:
        if |A[J1]| > AMAX:
            AMAX = |A[J1]|
            JMAX = IPIV + J
        J1 = J1 + LDA

    // ----- COMPUTE ROW MAXIMUM (RMAX) outside inner block -----
    // RMAX = max |A[IPIV, k]| for k in IEND_BLOCK+1..LIM
    // where LIM = NFRONT - KEEP(253) - NVSCHUR  (for PIVOT_OPTION=3)
    RMAX = 0
    for J = 1 to LIM - IEND_BLOCK:
        RMAX = max(RMAX, |A[IPIV row, column IEND_BLOCK+J]|)
        // Position: POSELT + (IEND_BLOCK+J-1)*LDA + (IPIV-1)

    // ----- NULL PIVOT CHECK -----
    if max(AMAX, RMAX, |PIVOT|) <= PIVNUL:
        // Entire row/column is essentially zero
        // Handle null pivot: set to FIXA or set row to zero with diagonal = 1
        goto accept_null_pivot

    // ----- 1x1 PIVOT TEST -----
    // After null check, extend RMAX to include abs(RMAX_NORELAX)
    RMAX = max(RMAX, |RMAX_NORELAX|)  // RMAX_NORELAX is 0 when PARPIV_T1=0
    if |PIVOT| >= UU * max(RMAX, AMAX) and |PIVOT| > max(SEUIL, tiny):
        NNEG++ if PIVOT < 0
        goto accept_1x1

    // ----- 2x2 PIVOT ATTEMPT -----
    // Need JMAX from the inner block search (AMAX > 0)
    if NPIVP1 == IEND_BLOCK: continue  // only 1 column, can't do 2x2
    if JMAX == 0: continue              // no off-diag found
    if max(|PIVOT|, RMAX, AMAX) <= tiny: continue

    // Refine RMAX: exclude the off-diagonal position (IPIV,JMAX)
    if RMAX < AMAX:
        // Recompute RMAX for column IPIV, excluding entry at JMAX
        for all entries in column IPIV outside inner block:
            RMAX = max(RMAX, |entry|)  // already done above
        // Also scan within-block entries excluding JMAX
        for JJ in [NPIVP1..IEND_BLOCK] row of col IPIV, JJ != JMAX:
            RMAX = max(RMAX, |A[JJ, IPIV]|)

    // ----- COMPUTE TMAX for candidate partner JMAX -----
    // TMAX = max |A[JMAX, k]| for all k != IPIV, k != JMAX
    APOSJ = POSELT + (JMAX-1)*LDA + NPIV
    POSPV2 = APOSJ + (JMAX - NPIVP1)   // diagonal of JMAX

    // Off-diagonal between IPIV and JMAX
    if IPIV < JMAX:
        OFFDAG = APOSJ + (IPIV - NPIVP1)    // A[IPIV, JMAX]
    else:
        OFFDAG = APOS + (JMAX - NPIVP1)     // A[JMAX, IPIV]

    TMAX = 0
    // Scan column JMAX in rows JMAX+1..LIM (outside), excluding IPIV
    // Scan column JMAX in rows NPIVP1..JMAX-1 (inside), excluding IPIV
    // (Detailed scanning code at lines 1553-1581)

    TMAX = max(TMAX, SEUIL/UU)  // minimum from TMAX_NORELAX

    // ----- 2x2 PIVOT TEST -----
    DETPIV = A[POSPV1] * A[POSPV2] - A[OFFDAG]^2

    // Reject if SEUIL check fails
    if SEUIL > 0 and sqrt(|DETPIV|) <= SEUIL: continue

    // Modified Bunch-Kaufman test for 2x2 block:
    //   |A[POSPV2]| * RMAX + AMAX * TMAX  <=  |DETPIV| / UU
    //   |A[POSPV1]| * TMAX + AMAX * RMAX  <=  |DETPIV| / UU
    if (|A[POSPV2]|*RMAX + AMAX*TMAX)*UU > |DETPIV|: continue
    if (|A[POSPV1]|*TMAX + AMAX*RMAX)*UU > |DETPIV|: continue

    // 2x2 pivot accepted!
    PIVSIZ = 2

    // Count negative eigenvalues:
    if DETPIV < 0:
        NNEG += 1       // one positive, one negative eigenvalue
    else if A[POSPV2] < 0:
        NNEG += 2       // both negative

    // Store DETPIV in sub-diagonal position
    A[POSELT + NPIV*(LDA+1) + 1] = DETPIV   // A[NPIV+2, NPIV+1] overwritten

    goto accept_2x2

// === END OF CANDIDATE LOOP (all candidates exhausted) ===
if IEND_BLOCK == NASS:
    INOPV = 1    // no pivot anywhere in NASS
else:
    INOPV = 2    // no pivot in this inner block, try extending
```

### Swap after pivot selection

After a pivot is found at position IPIV (and partner JMAX for 2x2):

```
accept:
    Inextpiv = max(NPIVP1 + PIVSIZ, IPIV + 1)

    // For 2x2: sort IPIV and JMAX
    for K = 1 to PIVSIZ:
        if PIVSIZ == 2:
            if K == 1: LPIV = min(IPIV, JMAX)
            if K == 2: LPIV = max(IPIV, JMAX)
        else:
            LPIV = IPIV

        if LPIV != NPIVP1:
            // Symmetric swap of rows/columns NPIVP1 and LPIV
            call DMUMPS_SWAP_LDLT(A, IW, IOLDPS, NPIVP1, LPIV, ...)

        NPIVP1 = NPIVP1 + 1

    if PIVSIZ == 2:
        // Store DETPIV at sub-diagonal of pivot block
        A[POSELT + NPIV*(LDA+1) + 1] = DETPIV
```

### The symmetric swap (DMUMPS_SWAP_LDLT)

Source: `dfac_front_aux.F` lines 2090-2151.

Swaps rows and columns `NPIVP1` and `IPIV` in a symmetric matrix. Note this is a SYMMETRIC swap — since `A` stores the full square, both the row and column must be swapped.

```
// Swap rows NPIVP1 and IPIV (in columns 1..NPIVP1-1)
// These are rows in the upper part, stored as columns:
dswap(NPIVP1-1, A[col NPIVP1, row 1], 1, A[col IPIV, row 1], 1)
// i.e., swap column NPIVP1 rows 1..NPIVP1-1 with column IPIV rows 1..NPIVP1-1

// Swap the "bridge" entries between NPIVP1 and IPIV
// Column NPIVP1+1 to IPIV-1 at row NPIVP1  with  column IPIV at rows NPIVP1+1..IPIV-1
dswap(IPIV - NPIVP1 - 1,
      A[POSELT + NPIVP1*LDA + (NPIVP1-1)], LDA,    // A[NPIVP1, NPIVP1+1], A[NPIVP1, NPIVP1+2], ...
      A[APOS + 1], 1)                               // A[NPIVP1+1, IPIV], A[NPIVP1+2, IPIV], ...

// Swap diagonals
swap(A[NPIVP1, NPIVP1],  A[IPIV, IPIV])

// Swap rows NPIVP1 and IPIV in columns IPIV+1..LASTROW2SWAP
if LASTROW2SWAP - IPIV > 0:
    dswap(LASTROW2SWAP - IPIV,
          A[APOS + LDA], LDA,       // A[NPIVP1, IPIV+1], stride LDA
          A[IDIAG + LDA], LDA)      // A[IPIV, IPIV+1], stride LDA

// Swap integer indices
swap(IW[row index NPIVP1], IW[row index IPIV])
swap(IW[col index NPIVP1], IW[col index IPIV])
```

Where `LASTROW2SWAP = NFRONT` (passed as `LIM_SWAP` in the caller, line 1642).

**Critical detail**: The swap operates on the FULL NFRONT rows/columns. For PIVOT_OPTION=3, LIM_SWAP = NFRONT, so every row beyond IPIV is also swapped.

---

## PHASE 2: WITHIN-PANEL UPDATE — FAC_MQ_LDLT

Source: `dfac_front_aux.F` lines 1677-1989.

Called immediately after each successful pivot. NPIV is the value BEFORE incrementing (the pivot just accepted is at position NPIV+1, or NPIV+1 and NPIV+2 for 2x2).

### Setup

```
NPIV_NEW = NPIV + PIVSIZ    // new pivot count after this step
NCB1 = LAST_ROW - IEND_BLOCK   // rows beyond the inner block to update
                                // For PIVOT_OPTION=3: LAST_ROW=NFRONT
NEL2 = IEND_BLOCK - NPIV_NEW   // remaining columns within this inner block

if NEL2 == 0:
    if IEND_BLOCK == NASS: IFINB = -1   // all of NASS done
    else: IFINB = 1                      // inner block done, panel continues
else:
    IFINB = 0                            // more pivots possible in this block
```

### 1x1 Pivot Update (PIVSIZ == 1)

The pivot is at diagonal position `(NPIV+1, NPIV+1)`.

```
APOS = POSELT + NPIV*(LDA + 1)    // A[NPIV+1, NPIV+1] (0-based offset)
VALPIV = 1.0 / A[APOS]            // reciprocal of pivot

// Process rows NPIV+2 through NPIV+2+NEL2+NCB1-1
// i.e., all rows from NPIV+2 to LAST_ROW

// For each row I = 1 to NEL2 (within inner block, triangular part):
for I = 1 to NEL2:
    K1POS = APOS + I * LDA    // position A[NPIV+1, NPIV+1+I] = A[pivot_row, col NPIV+1+I]
                               // Actually this is column NPIV+1+I at row NPIV+1

    // STEP 1: Save original L entry to U storage (pivot row)
    A[APOS + I] = A[K1POS]    // Copy to pivot row: A[NPIV+1+I, NPIV+1] -> A[NPIV+1, NPIV+1+I]
                               // (This stores the UNSCALED original value in the upper triangle)

    // STEP 2: Scale L entry by 1/pivot
    A[K1POS] = A[K1POS] * VALPIV    // L[NPIV+1+I, NPIV+1] = A[...] / d

    // STEP 3: Update trailing entries (triangular within block)
    // Only update columns 1..I (triangular: j <= I)
    for JJ = 1 to I:
        A[K1POS + JJ] = A[K1POS + JJ] - A[K1POS] * A[APOS + JJ]
        // A[NPIV+1+JJ, NPIV+1+I] -= L[NPIV+1+I, NPIV+1] * U[NPIV+1, NPIV+1+JJ]
        // where U[NPIV+1, NPIV+1+JJ] = A[APOS+JJ] = original (unscaled) value saved in step 1
    endfor
endfor

// For each row I = NEL2+1 to NEL2+NCB1 (beyond inner block, rectangular part):
for I = NEL2+1 to NEL2+NCB1:
    K1POS = APOS + I * LDA

    // STEP 1: Save original to U
    A[APOS + I] = A[K1POS]

    // STEP 2: Scale L entry
    A[K1POS] = A[K1POS] * VALPIV

    // STEP 3: Rectangular update (all NEL2 columns)
    for JJ = 1 to NEL2:
        A[K1POS + JJ] = A[K1POS + JJ] - A[K1POS] * A[APOS + JJ]
    endfor
endfor
```

**Storage layout after 1x1 update**: The pivot column of the frontal matrix now has:
- `A[APOS]` = original `d` value (the 1x1 pivot, NOT its reciprocal)
- `A[APOS + I]` for I=1..NEL2+NCB1 = **original unscaled** sub-diagonal values (the "U" copy)
- Column NPIV+1, rows NPIV+2..LAST_ROW: the **scaled** L entries: `L[i, NPIV+1] = a[i, NPIV+1] / d`

The Schur complement update formula is:
```
A[row, col] -= L[row, pivotcol] * U[pivotcol, col]
             = (a_orig / d) * a_orig_col
```

This is equivalent to `A -= L * D * L^T` because:
- `L[i,k] = a[i,k] / d[k]`
- `U_saved[k,j] = a[k,j]` (the original, which equals `D[k]*L[j,k]` since `L[j,k] = a[j,k]/d[k]`)
- So the update is `A[i,j] -= L[i,k] * D[k] * L[j,k] = (a[i,k]/d[k]) * a[j,k]`

### 2x2 Pivot Update (PIVSIZ == 2)

The 2x2 pivot block occupies positions `(NPIV+1, NPIV+1)`, `(NPIV+2, NPIV+2)`, and the off-diagonal at `(NPIV+2, NPIV+1)`.

```
POSPV1 = POSELT + NPIV*(LDA + 1)         // A[NPIV+1, NPIV+1]
POSPV2 = POSPV1 + LDA + 1                // A[NPIV+2, NPIV+2]
OFFDAG_OLD = POSPV2 - 1                  // A[NPIV+1, NPIV+2] = A[row NPIV+1, col NPIV+2]
OFFDAG = POSPV1 + 1                      // A[NPIV+2, NPIV+1] = A[row NPIV+2, col NPIV+1]

// Read pivot block values
SWOP = A[POSPV2]         // d22 (diagonal of NPIV+2)
DETPIV = A[OFFDAG]        // This was stored by FAC_I_LDLT: det(D) = d11*d22 - d12^2

// Compute D^{-1} elements (scaled by 1/det):
A22 = A[POSPV1] / DETPIV    // d11 / det  (this gives (D^{-1})_22)
A11 = SWOP / DETPIV          // d22 / det  (this gives (D^{-1})_11)
A12 = -A[OFFDAG_OLD] / DETPIV // -d12 / det (this gives (D^{-1})_12)

// Fix storage: move off-diagonal from upper to lower position
A[OFFDAG] = A[OFFDAG_OLD]    // A[NPIV+2, NPIV+1] = d12 (the off-diagonal)
A[OFFDAG_OLD] = 0             // A[NPIV+1, NPIV+2] = 0 (cleared)

// Starting positions for L columns
LPOS1 = POSPV2 + LDA - 1   // A[NPIV+1, NPIV+3] = row NPIV+1 of column NPIV+3
LPOS2 = LPOS1 + 1           // A[NPIV+2, NPIV+3]

// === NON-VE (standard CPU) PATH ===
// Process within-block rows (triangular update): J2 = 1..NEL2
JJ = POSPV2 + LDA - 1    // A[NPIV+1, NPIV+3]

for J2 = 1 to NEL2:
    K1 = JJ              // A[NPIV+1, column NPIV+2+J2]
    K2 = JJ + 1          // A[NPIV+2, column NPIV+2+J2]

    // Compute L * D^{-1} entries (negated for subtraction)
    MULT1 = -(A11*A[K1] + A12*A[K2])
    MULT2 = -(A12*A[K1] + A22*A[K2])

    // Save original sub-diagonal values to U storage (pivot rows)
    A[POSPV1 + 2 + (J2-1)] = A[K1]   // U row 1: original a[NPIV+1, col]
    A[POSPV2 + 1 + (J2-1)] = A[K2]   // U row 2: original a[NPIV+2, col]

    // Update trailing (triangular): columns NPIV+3..NPIV+2+J2
    K1 = POSPV1 + 2       // start of saved U row 1
    K2 = POSPV2 + 1       // start of saved U row 2
    for IROW = IBEG to IEND:   // triangular region
        A[IROW] = A[IROW] + MULT1*A[K1] + MULT2*A[K2]
        K1++; K2++
    endfor

    // Store scaled L entries (L*D^{-1}) back in the L position
    A[JJ]     = -MULT1    // L[NPIV+1, col] * D^{-1} part 1
    A[JJ + 1] = -MULT2    // L[NPIV+2, col] * D^{-1} part 2

    JJ += LDA
    // IBEG and IEND advance to cover the growing triangular region
endfor

// Process beyond-block rows (rectangular update): rows IEND_BLOCK+1..LAST_ROW
for J2 = 1 to LAST_ROW - IEND_BLOCK:
    K1 = JJ_LOC           // A[NPIV+1, column IEND_BLOCK+J2]
    K2 = JJ_LOC + 1       // A[NPIV+2, column IEND_BLOCK+J2]

    MULT1 = -(A11*A[K1] + A12*A[K2])
    MULT2 = -(A12*A[K1] + A22*A[K2])

    // Save originals to U
    A[POSPV1 + 2 + NEL2 + (J2-1)] = A[K1]
    A[POSPV2 + 1 + NEL2 + (J2-1)] = A[K2]

    // Rectangular update (all NEL2 columns)
    K1 = POSPV1 + 2; K2 = POSPV2 + 1
    for IROW = IBEG_LOC to IEND_LOC:
        A[IROW] = A[IROW] + MULT1*A[K1] + MULT2*A[K2]
        K1++; K2++
    endfor

    // Store scaled L
    A[JJ_LOC]     = -MULT1
    A[JJ_LOC + 1] = -MULT2
endfor
```

**Critical detail about the 2x2 update formula**: The L entries stored back are `L * D^{-1}`, NOT simply `L / d` as in the 1x1 case. Specifically:

```
[L1]     [A11  A12] [a1]     [D^{-1}] [original_col]
[L2]  =  [A12  A22] [a2]  =                          
```

where `(A11, A12, A22)` are the elements of `D^{-1}` and `(a1, a2)` are the original off-diagonal entries. The stored L entries are `D^{-1} * original`, and the update uses `-L * original^T`.

**Storage after 2x2 update**:
- `A[POSPV1]` = original `d11` (diagonal of first pivot)
- `A[POSPV2]` = original `d22` (diagonal of second pivot)
- `A[OFFDAG]` = `A[NPIV+2, NPIV+1]` = original `d12` (off-diagonal)
- `A[OFFDAG_OLD]` = `A[NPIV+1, NPIV+2]` = 0 (zeroed)
- Pivot row 1 (`A[POSPV1+2..POSPV1+2+NEL2+NCB1-1]`): original unscaled `a[NPIV+1, col]` values
- Pivot row 2 (`A[POSPV2+1..POSPV2+1+NEL2+NCB1-1]`): original unscaled `a[NPIV+2, col]` values
- L column 1 (col NPIV+3 onwards, row NPIV+1): `D^{-1}_{11}*a1 + D^{-1}_{12}*a2`
- L column 2 (col NPIV+3 onwards, row NPIV+2): `D^{-1}_{12}*a1 + D^{-1}_{22}*a2`

### IFINB return values

```
IFINB = 0   -> more columns remain in this inner block, continue pivot-by-pivot
IFINB = 1   -> inner block exhausted (IEND_BLOCK reached), but more in BLR panel
IFINB = -1  -> NASS reached (all fully-summed columns eliminated or delayed)
```

---

## PHASE 3: BETWEEN-PANEL BLAS-3 UPDATE — FAC_SQ_LDLT

Source: `dfac_front_aux.F` lines 1990-2089.

Called in two places:
1. After each inner block completes (within a BLR panel): updates columns IEND_BLOCK+1..IEND_BLR
2. After each BLR panel completes: updates columns IEND_BLR+1..NASS and rows NASS+1..NFRONT

### What it computes

This routine applies accumulated pivots from `IBEG_BLOCK..NPIV` (the block just factored) to trailing columns. It has two BLAS-3 operations:

**TRSM**: Solve `U_block * X = A_target` where U_block is the upper triangular factor from the pivot block. This computes the L entries for the trailing columns.

**GEMM**: Update the Schur complement using the newly computed L entries.

### Detailed algorithm

```
NPIV_BLOCK = NPIV - IBEG_BLOCK + 1   // number of pivots in this block
if NPIV_BLOCK == 0: return            // nothing to do

// Parameters (for the within-panel call from the driver):
//   IBEG_BLOCK, IEND_BLOCK: the inner block range
//   NPIV: current pivot count
//   FIRST_ROW_TRSM = -6666, LAST_ROW_TRSM = -6666 (TRSM not called)
//   LAST_COL_GEMM = IEND_BLR
//   LAST_ROW_GEMM = LAST_ROW (= NFRONT for PIVOT_OPTION=3)
//   CALL_TRSM = false, CALL_GEMM = true

// For the inter-panel call:
//   FIRST_ROW_TRSM = IEND_BLR
//   LAST_ROW_TRSM = NASS
//   LAST_COL_GEMM = NASS
//   LAST_ROW_GEMM = LAST_ROW (= NFRONT)
//   CALL_TRSM = (PIVOT_OPTION <= 1), CALL_GEMM = true

NEL1 = LAST_COL_GEMM - IEND_BLOCK    // number of columns to update

// === TRSM (if CALL_TRSM) ===
if CALL_TRSM and NEL1 != 0:
    APOS = POSELT + (IBEG_BLOCK-1)*LDA + (IBEG_BLOCK-1)  // pivot block diagonal
    LPOS = POSELT + FIRST_ROW_TRSM*LDA + (IBEG_BLOCK-1)  // L target
    UPOS = POSELT + (IBEG_BLOCK-1)*LDA + FIRST_ROW_TRSM  // U target

    NRHS_TRSM = LAST_ROW_TRSM - FIRST_ROW_TRSM

    // Triangular solve: L * X = B  where L is the upper-triangular
    // pivot block (unit diagonal, stored transposed)
    dtrsm('L', 'U', 'T', 'U', NPIV_BLOCK, NRHS_TRSM,
          ONE, A[APOS], LDA, A[LPOS], LDA)

    // Copy L to U and scale by D^{-1}
    call DMUMPS_FAC_LDLT_COPY2U_SCALEL(NRHS_TRSM, 1, KEEP(424),
         NFRONT, NPIV_BLOCK, ..., LPOS, UPOS, APOS)

// === GEMM ===
if CALL_GEMM and NEL1 != 0:
    // Part 1: Symmetric update of [IEND_BLOCK+1..LAST_COL_GEMM] x [IEND_BLOCK+1..LAST_COL_GEMM]
    // This is a symmetric rank-k update, done with blocked DGEMM
    for IROW = IEND_BLOCK+1 to LAST_COL_GEMM, step BLSIZE:
        Block = min(BLSIZE, LAST_COL_GEMM - IROW + 1)
        LPOS = POSELT + (IROW-1)*LDA + (IBEG_BLOCK-1)
        UPOS = POSELT + (IBEG_BLOCK-1)*LDA + (IROW-1)
        APOS = POSELT + (IROW-1)*LDA + (IROW-1)

        dgemm('N', 'N', Block, LAST_COL_GEMM - IROW + 1, NPIV_BLOCK,
              -1.0, A[UPOS], LDA, A[LPOS], LDA, 1.0, A[APOS], LDA)
    endfor

    // Part 2: Rectangular update of CB rows
    if LAST_ROW_GEMM > LAST_COL_GEMM:
        LPOS = POSELT + LAST_COL_GEMM*LDA + (IBEG_BLOCK-1)
        UPOS = POSELT + (IBEG_BLOCK-1)*LDA + IEND_BLOCK
        APOS = POSELT + LAST_COL_GEMM*LDA + IEND_BLOCK

        dgemm('N', 'N', NEL1, LAST_ROW_GEMM - LAST_COL_GEMM, NPIV_BLOCK,
              -1.0, A[UPOS], LDA, A[LPOS], LDA, 1.0, A[APOS], LDA)
```

### The COPY2U_SCALEL operation

Source: lines 2152-2265.

This is called from FAC_SQ_LDLT after TRSM. It processes each pivot column:

```
for I = 1 to NPIV_BLOCK:
    DPOS = A_DPOS + (I-1)*LDA + (I-1)   // diagonal of pivot I

    // Check if this is a 2x2 pivot (by checking sign of IW column index)
    PIVOT_2X2 = (IW[OFFSET_IW + I - 1] <= 0)
    // But skip the SECOND column of a 2x2 (handled together with first)
    if not PIVOT_2X2 and I > 1:
        if IW[OFFSET_IW + I - 2] <= 0: cycle  // this is col 2 of a 2x2

    if not PIVOT_2X2:
        // 1x1: copy L to U, then scale L by 1/d
        A11 = 1.0 / A[DPOS]
        for J = 1 to Block2:  // Block2 = number of trailing rows to process
            A[UPOS + ...] = A[LPOS + ...]    // copy to U (transposed position)
        endfor
        for J = 1 to Block2:
            A[LPOS + ...] = A[LPOS + ...] * A11   // scale L by 1/d
        endfor
    else:
        // 2x2: copy both columns to U, then compute L*D^{-1}
        dcopy(Block2, L_col1, LDA, U_row1, 1)
        dcopy(Block2, L_col2, LDA, U_row2, 1)

        A11 = A[POSPV1]; A22 = A[POSPV2]; A12 = A[OFFDAG]
        DETPIV = A11*A22 - A12^2
        // Compute D^{-1}:
        A22_inv = A11 / DETPIV
        A11_inv = A[POSPV2] / DETPIV
        A12_inv = -A12 / DETPIV

        for J = 1 to Block2:
            MULT1 = A11_inv * L[J, col1] + A12_inv * L[J, col2]
            MULT2 = A12_inv * L[J, col1] + A22_inv * L[J, col2]
            L[J, col1] = MULT1    // overwrite L with L*D^{-1}
            L[J, col2] = MULT2
        endfor
```

**Critical insight**: After TRSM + COPY2U_SCALEL:
- The **U** region (upper triangle, i.e., pivot rows) contains the **original unscaled** values
- The **L** region (below pivot rows) contains `L * D^{-1}` (L entries pre-multiplied by D inverse)
- The GEMM then computes `A -= U^T * L` which equals `A -= (D*L^T) * (L*D^{-1}) * ... ` wait, no.

Let me re-check. The GEMM call is:
```
dgemm('N', 'N', ..., -1.0, A[UPOS], LDA, A[LPOS], LDA, 1.0, A[APOS], LDA)
```

So it computes `A_trailing -= U * L` where:
- `U[k, j]` = original `a[k, j]` (the unscaled saved values) — these are `D * L^T` entries
- `L[i, k]` = `L[i,k] * D^{-1}[k,k]` (the scaled L entries)

Wait, that gives `A -= (D*L^T)^T * (L*D^{-1})` ... No. Let me reconsider.

The positions: UPOS points to pivot rows (the "U" copy), LPOS points to L entries below the panel. The GEMM multiplies columns of U (at UPOS) with rows of L (at LPOS).

Actually re-reading the GEMM call more carefully: `dgemm('N','N', Block, M, NPIV_BLOCK, -1.0, A(UPOS), LDA, A(LPOS), LDA, 1.0, A(APOS), LDA)`. 

UPOS = `POSELT + (IBEG_BLOCK-1)*LDA + (IROW-1)` — this is the U storage, which is in the pivot rows, at column positions of the target. The layout is: for the k-th pivot (at row IBEG_BLOCK-1+k), the entries in columns IROW..LAST are stored. These contain the original unscaled values.

LPOS = `POSELT + (IROW-1)*LDA + (IBEG_BLOCK-1)` — this is below the pivot block, at the target rows. These contain `L * D^{-1}`.

So `A -= U_block * L_block` where:
- `U_block` dimensions: NPIV_BLOCK x (LAST_COL_GEMM - IROW + 1) — these are `D * L^T` (original values)
- `L_block` dimensions: NPIV_BLOCK x Block — these are `L * D^{-1}`

The update: `A[i,j] -= sum_k U[k,j] * L[k,i]` ... Hmm, with `dgemm('N','N')`, the result is `C = alpha * A * B + beta * C`, so:

`A_trailing -= U_block^{at UPOS} * L_block^{at LPOS}`

where U_block is `Block x NPIV_BLOCK` (UPOS has leading dimension LDA, first dim = Block = columns of target in the trailing region) and L_block is `NPIV_BLOCK x (LAST_COL-IROW+1)`.

Wait, I need to re-check the exact argument order:
```
dgemm('N', 'N', Block, LAST_COL_GEMM - IROW + 1, NPIV_BLOCK,
      -1.0, A(UPOS), LDA, A(LPOS), LDA, 1.0, A(APOS), LDA)
```

- M = Block, N = LAST_COL_GEMM - IROW + 1, K = NPIV_BLOCK
- C(M,N) = alpha * A(M,K) * B(K,N) + beta * C(M,N)
- A at UPOS: M=Block rows x K=NPIV_BLOCK cols
- B at LPOS: K=NPIV_BLOCK rows x N=(LAST_COL-IROW+1) cols

UPOS = `POSELT + (IBEG_BLOCK-1)*LDA + (IROW-1)` — row IROW, column IBEG_BLOCK in the matrix. So the A matrix of DGEMM has rows=IROW..IROW+Block-1, columns=IBEG_BLOCK..IBEG_BLOCK+NPIV_BLOCK-1. These are the **U entries** (original values stored in the upper triangle / pivot rows).

LPOS = `POSELT + (IROW-1)*LDA + (IBEG_BLOCK-1)` — column IROW, row IBEG_BLOCK. So the B matrix has rows=IBEG_BLOCK..IBEG_BLOCK+NPIV_BLOCK-1, columns=IROW..IROW+N-1. These are the **L entries** (scaled by D^{-1}).

So: `C[i,j] -= U[i,k] * L[k,j]` where i is in target rows IROW..IROW+Block-1, j is in target columns IROW..LAST_COL, and k ranges over IBEG_BLOCK..NPIV.

This means U contains the entries at `A[target_row, pivot_col]` (the row of the original unscaled values) and L contains `A[pivot_row, target_col]` which has been scaled by `D^{-1}`.

The mathematical update is:
```
A_trailing[i,j] -= sum_k  original[i,k] * (D^{-1} * original)[k,j]
                 = sum_k  (L^T * D)[i,k] * (D^{-1} * L^T)^T ... 
```

This is the standard `A -= L_orig * D^{-1} * L_orig^T` Schur complement update, split as `(L_orig) * (D^{-1} * L_orig^T)` where:
- U stores `L_orig^T` (the original unscaled values, transposed into the upper triangle)
- L stores `D^{-1} * L_orig` (but stored column-major in the lower part, so it's effectively `L_orig * D^{-1}` when read by column)

The GEMM `C -= A * B` with `A = U_entries` and `B = L_entries` gives:

`A_trailing[i,j] -= U[i,:] * L[:,j] = sum_k original[k,i]^T * (D^{-1} * original)^T[k,j]`

No wait. Let me just be very precise about what's stored where. After COPY2U_SCALEL for a 1x1 pivot `d` at position `p`:
- U row: `A[p, j]` for `j > p` = `original_a[p, j]` = the value that was at position `(p, j)` before this pivot
- L col: `A[i, p]` for `i > p` = `original_a[i, p] / d` = the L factor

For the GEMM, U entries at UPOS are rows of the target blocks intersected with pivot columns. The A matrix argument reads `A[target_row, pivot_col]` = `original[target_row, pivot_col]`.

The B matrix reads `A[pivot_row, target_col]` = `original[pivot_row, target_col] / d` = `L[target_col, pivot_row]^T / d`... Hmm, actually for the LDLT the L entries are in the COLUMNS below the diagonal, and by symmetry `A[i,k]` (row i, col k with i > k) is the same as `A[k,i]`. After the COPY2U step, the column entries `A[i, pivot]` are scaled by `1/d`, and the row entries `A[pivot, j]` have the original unscaled values.

So the GEMM computes:
```
trailing[i,j] -= sum_k  A[i, pivot_k] * A[pivot_k, j]
              -= sum_k  (original[i,k] / d_k) * original[j,k]   ... wait, A[pivot_k, j] is in U
```

Hmm no. A[pivot_k, j] is at UPOS which stores the original value. And A[i, pivot_k] is at LPOS which stores the value divided by d. But these are accessed by column, so:

Actually, looking at the DGEMM arguments again:
- First matrix (M x K) is at UPOS = `A[IROW, IBEG_BLOCK]` — these are the **saved originals** in the U region (row entries of the pivot rows stored into what was the upper triangle).
- Second matrix (K x N) is at LPOS = `A[IBEG_BLOCK, IROW]` — these are the **scaled L** entries (column entries below the diagonal, scaled by D^{-1}).

But wait — UPOS points to `row = IROW` through `row = IROW + Block - 1` in `columns = IBEG_BLOCK` through `IBEG_BLOCK + NPIV_BLOCK - 1`. By symmetry in the original matrix, `A[IROW, IBEG_BLOCK]` is an entry in the lower triangle (IROW > IBEG_BLOCK for the trailing region). But after COPY2U_SCALEL, this position now holds the original unscaled copy (U storage).

And LPOS points to `column = IROW` through `column = IROW + N - 1` in `rows = IBEG_BLOCK` through `IBEG_BLOCK + NPIV_BLOCK - 1`. These are the pivot rows at the trailing columns. After COPY2U_SCALEL, the L entries are at column positions (below pivot rows), but wait — `A[IBEG_BLOCK, IROW]` is in the upper triangle (IBEG_BLOCK < IROW). In the factored region, this is the L column.

Actually, I think I was overcomplicating this. The key insight from the MQ routine is clearer: during the rank-1 update,

For 1x1: `A[APOS+I]` stores the original value (call it `u_i`), and `A[K1POS]` after scaling stores `l_i = u_i / d`. The update is `A[row,col] -= l_row * u_col`. This is mathematically:

```
A -= l * u^T = (1/d * a_col) * a_col^T
```

Which is exactly `A -= L[:,k] * D[k,k] * L[:,k]^T` because `u = D * L^T`.

So in the GEMM: `C -= A_upos * B_lpos` computes `trailing -= original * scaled` which is `trailing -= (D * L^T) * L` = `L^T * D * L`. This is correct for the Schur complement.

---

## PHASE 4: COMPLETE PSEUDOCODE

Here is the full implementable algorithm, cleaned up and using 0-based indexing for the Rust implementation.

### Notation for 0-based

- Matrix M is NFRONT x NFRONT, column-major, stored as `M[row + col * LDA]`
- Pivot count starts at 0

```rust
// === INITIALIZATION ===
let mut npiv: usize = 0;  // number of eliminated pivots
let nbkjib = inner_block_size(nass);  // typically 16 or 32
let nblr = 4 * 32;  // KEEP(420) = 128

let mut iend_blr: usize = 0;
let mut iend_block: usize = 0;
let mut last_panel = false;

// === OUTER PANEL LOOP ===
while iend_blr < nass && !last_panel {
    iend_blr = min(iend_blr + nblr, nass);

    // === INNER BLOCK LOOP ===
    while iend_block < iend_blr && !last_panel {
        let ibeg_block = npiv;  // 0-based start of this inner block
        iend_block = min(iend_block + nbkjib, iend_blr);

        // === PIVOT LOOP ===
        loop {
            // STEP 1: PIVOT SEARCH (FAC_I_LDLT)
            let (inopv, pivsiz, swap_info) = find_pivot_ldlt(
                &mut A, nfront, nass, npiv,
                ibeg_block, iend_block,  // 0-based range [ibeg_block, iend_block)
                uu, seuil,
            );

            match inopv {
                1 => {  // No pivot found in entire NASS range
                    // Try static pivoting if enabled
                    if static_mode {
                        static_pivot(&mut A, npiv, nfront, cseuil);
                        // pivsiz = 1, continue
                    } else {
                        last_panel = true;
                        break;
                    }
                }
                2 => break,  // No pivot in this inner block, advance
                _ => {}       // Pivot found, continue
            }

            // STEP 2: SWAP PIVOT INTO POSITION
            // (Already done inside find_pivot_ldlt via SWAP_LDLT)

            // STEP 3: WITHIN-PANEL UPDATE (FAC_MQ_LDLT)
            let last_row = nfront;  // PIVOT_OPTION=3
            let nvschur_k253 = nvschur + keep_253;
            let ifinb = update_within_block_ldlt(
                &mut A, nfront, nass, npiv, pivsiz,
                iend_block, last_row, nvschur_k253,
            );

            // STEP 4: Mark 2x2 pivot in metadata
            if pivsiz == 2 {
                // Negate column index of second pivot position
                col_indices[npiv + 1] = -col_indices[npiv + 1];
            }

            npiv += pivsiz;

            match ifinb {
                0  => continue,         // more pivots in this block
                1  => break,            // block done, panel continues
                -1 => { last_panel = true; break; }  // NASS done
            }
        }

        // STEP 5: BETWEEN-BLOCK BLAS-3 UPDATE (within BLR panel)
        if iend_blr > iend_block {
            let npiv_block = npiv - ibeg_block;
            if npiv_block > 0 {
                // No TRSM here (handled by MQ for PIVOT_OPTION=3)
                // GEMM: update columns iend_block..iend_blr-1 and rows iend_block..nfront-1
                between_block_update_ldlt(
                    &mut A, nfront, ibeg_block, iend_block, npiv,
                    iend_blr,   // LAST_COL_GEMM (0-based: columns iend_block..iend_blr-1)
                    nfront,     // LAST_ROW_GEMM
                );
            }
        }
    }  // end inner block loop

    // STEP 6: INTER-PANEL BLAS-3 UPDATE
    let ibeg_blr = ... ;  // start of this panel
    let npiv_block = npiv - ibeg_blr;
    if npiv_block > 0 {
        // TRSM + COPY2U_SCALEL + GEMM for columns nass..nfront-1
        inter_panel_update_ldlt(
            &mut A, nfront, nass,
            ibeg_blr, iend_blr, npiv,
            nfront,  // LAST_ROW
        );
    }
}
```

### find_pivot_ldlt (corresponds to FAC_I_LDLT)

```rust
fn find_pivot_ldlt(
    A: &mut [f64], lda: usize, nass: usize, npiv: usize,
    ibeg_block: usize, iend_block: usize,  // 0-based half-open range
    uu: f64, seuil: f64,
) -> (i32, usize, Option<SwapInfo>) {
    // inopv: 0=found, 1=none in NASS, 2=none in this block
    // pivsiz: 1 or 2
    let pivnul = dkeep_1;  // DKEEP(1), null pivot threshold

    for ipiv in ibeg_block..iend_block {  // try each candidate column
        let pivot = A[ipiv + ipiv * lda];  // diagonal

        // SCAN AMAX, JMAX: max off-diagonal in column IPIV within [npiv, iend_block)
        let (amax, jmax) = scan_column_within_block(A, lda, npiv, ipiv, iend_block);

        // SCAN RMAX: max |entry| in column IPIV outside [npiv, iend_block)
        // For PIVOT_OPTION=3: rows iend_block..nfront-1 (excluding NVSCHUR/K253)
        let rmax = scan_column_outside_block(A, lda, ipiv, iend_block, lim);

        // NULL PIVOT CHECK
        if max3(amax, rmax, pivot.abs()) <= pivnul {
            handle_null_pivot(A, lda, ipiv, npiv, fixa);
            return (0, 1, ...);
        }

        // 1x1 TEST: |d_ii| >= uu * max(rmax, amax)  AND  |d_ii| > seuil
        if pivot.abs() >= uu * f64::max(rmax, amax)
           && pivot.abs() > f64::max(seuil, f64::MIN_POSITIVE) {
            if pivot < 0.0 { nneg += 1; }
            // SWAP ipiv into position npiv
            if ipiv != npiv {
                symmetric_swap(A, lda, nfront, npiv, ipiv);
            }
            return (0, 1, swap_info);
        }

        // 2x2 ATTEMPT
        if jmax.is_none() || npiv + 1 == iend_block { continue; }
        let jmax = jmax.unwrap();

        // Refine RMAX for column IPIV (exclude JMAX entry)
        // ... (re-scan excluding the off-diagonal partner)

        // SCAN TMAX: max |entry| in column JMAX, excluding IPIV entry
        let tmax = scan_column_for_partner(A, lda, npiv, jmax, ipiv, lim);
        let tmax = f64::max(tmax, seuil / uu);

        let d_ii = A[ipiv + ipiv * lda];
        let d_jj = A[jmax + jmax * lda];
        let d_ij = if ipiv < jmax {
            A[ipiv + jmax * lda]
        } else {
            A[jmax + ipiv * lda]
        };

        let detpiv = d_ii * d_jj - d_ij * d_ij;
        let abs_det = detpiv.abs();

        // SEUIL check on 2x2
        if seuil > 0.0 && abs_det.sqrt() <= seuil { continue; }

        // MODIFIED BUNCH-KAUFMAN TEST:
        // Both conditions must hold:
        //   (|d_jj| * rmax + amax * tmax) * uu <= |det|
        //   (|d_ii| * tmax + amax * rmax) * uu <= |det|
        if (d_jj.abs() * rmax + amax * tmax) * uu > abs_det { continue; }
        if (d_ii.abs() * tmax + amax * rmax) * uu > abs_det { continue; }

        // 2x2 ACCEPTED
        // Count negative eigenvalues
        if detpiv < 0.0 {
            nneg += 1;   // one positive + one negative
        } else if d_jj < 0.0 {
            nneg += 2;   // both negative
        }

        // SWAP: put min(ipiv,jmax) at npiv, max(ipiv,jmax) at npiv+1
        let first = min(ipiv, jmax);
        let second = max(ipiv, jmax);
        if first != npiv {
            symmetric_swap(A, lda, nfront, npiv, first);
        }
        if second != npiv + 1 {
            symmetric_swap(A, lda, nfront, npiv + 1, second);
        }

        // Store DETPIV at sub-diagonal of pivot block
        A[npiv + 1 + npiv * lda] = detpiv;  // A[npiv+1, npiv] in 0-based

        return (0, 2, swap_info);
    }

    // No pivot found
    if iend_block == nass {
        return (1, 0, None);  // exhausted all of NASS
    } else {
        return (2, 0, None);  // only this inner block exhausted
    }
}
```

### symmetric_swap (corresponds to DMUMPS_SWAP_LDLT)

```rust
fn symmetric_swap(A: &mut [f64], lda: usize, nfront: usize, p: usize, q: usize) {
    // Swap rows and columns p and q in the symmetric matrix (0-based)
    // p < q always
    assert!(p < q);

    // 1. Swap entries in columns 0..p-1 (rows p and q)
    //    A[p, 0..p-1]  <-> A[q, 0..p-1]
    for col in 0..p {
        A.swap(p + col * lda, q + col * lda);
    }

    // 2. Swap the "bridge" region:
    //    A[p, p+1..q-1]  <-> A[p+1..q-1, q]
    //    These are: row p in columns p+1..q-1 vs column q in rows p+1..q-1
    for k in 1..(q - p) {
        // A[p, p+k] (stored at row p, column p+k) is entry (p, p+k) in upper triangle
        // But since symmetric, the lower entry A[p+k, p] is the one with stride LDA
        // MUMPS stores FULL matrix, so:
        // Swap A[p, p+k] (= A[p + (p+k)*lda]) with A[p+k, q] (= A[p+k + q*lda])
        // Wait - re-reading MUMPS code:
        //   dswap(IPIV-NPIVP1-1,
        //         A(POSELT+NPIVP1*LDA8 + NPIVP1-1), LDA,   // col NPIVP1+1 row NPIVP1, stride LDA
        //         A(APOS + 1), 1)                            // col q rows p+1..q-1, stride 1
        // 1-based: A[NPIVP1, NPIVP1+1] stride LDA vs A[NPIVP1+1, IPIV] stride 1
        // The first is row p (0-based) in columns p+1, p+2, ..., q-1
        //   accessed as A[p + (p+1)*lda], A[p + (p+2)*lda], ... stride = lda
        // The second is column q (0-based), rows p+1, p+2, ..., q-1
        //   accessed as A[(p+1) + q*lda], A[(p+2) + q*lda], ... stride = 1
        let idx1 = p + (p + k) * lda;
        let idx2 = (p + k) + q * lda;
        A.swap(idx1, idx2);
    }

    // 3. Swap diagonals
    A.swap(p + p * lda, q + q * lda);

    // 4. Swap entries in columns q+1..nfront-1 (rows p and q)
    //    A[p, q+1..nfront-1] <-> A[q, q+1..nfront-1]
    for col in (q + 1)..nfront {
        // row p in column col vs row q in column col
        let idx1 = p + col * lda;
        let idx2 = q + col * lda;
        A.swap(idx1, idx2);
    }

    // 5. Swap row/column index metadata
    // swap(row_indices[p], row_indices[q])
    // swap(col_indices[p], col_indices[q])
}
```

### update_within_block_ldlt (corresponds to FAC_MQ_LDLT)

```rust
fn update_within_block_ldlt(
    A: &mut [f64], lda: usize, nass: usize, npiv: usize, pivsiz: usize,
    iend_block: usize, last_row: usize, nvschur_k253: usize,
) -> i32 {
    let npiv_new = npiv + pivsiz;
    let nel2 = iend_block - npiv_new;       // remaining cols in inner block
    let ncb1 = last_row - iend_block;       // rows beyond inner block

    let ifinb = if nel2 == 0 {
        if iend_block == nass { -1 } else { 1 }
    } else {
        0
    };

    if pivsiz == 1 {
        // === 1x1 PIVOT UPDATE ===
        let apos = npiv + npiv * lda;        // diagonal A[npiv, npiv]
        let d = A[apos];
        let valpiv = 1.0 / d;

        // Process each trailing row
        // "I" ranges from 1 to nel2 + ncb1
        // Row index (0-based) = npiv + I = npiv+1, npiv+2, ...
        // The entry at A[row, npiv] is at: row + npiv * lda
        // But in column-major with the transposed view used by MUMPS:
        // A[npiv, col] for col > npiv is in column col at row npiv
        //
        // MUMPS convention: APOS = A[npiv, npiv] (0-based)
        // LPOS = APOS + I * LDA = A[npiv, npiv + I] = A[row npiv, col npiv+I]
        //
        // In the symmetric matrix, A[npiv, npiv+I] = A[npiv+I, npiv]
        // These are accessed via columns: column npiv+I at row npiv
        //
        // Wait: rechecking. APOS = POSELT + NPIV*(NFRONT+1) in 1-based.
        // That's row=NPIV+1, col=NPIV+1 in 1-based = row=npiv, col=npiv in 0-based.
        // LPOS = APOS + I * LDA = row npiv, col npiv+I
        //
        // A[K1POS] = A[npiv, npiv+I] in the FULL matrix
        // This is the (npiv+I)-th column at row npiv = the off-diagonal
        // By symmetry this should equal A[npiv+I, npiv]
        // After the update: A[K1POS] = L[npiv+I, npiv] = old / d

        for i in 1..=(nel2 + ncb1) {
            let col = npiv + i;       // 0-based column index
            let k1pos = apos + i * lda;  // A[npiv, col] = A[npiv, npiv+i]

            // Save original to "U" storage (pivot row)
            A[apos + i] = A[k1pos];   // A[npiv+i, npiv] gets the original

            // Scale: L[col, npiv] = original / d
            A[k1pos] = A[k1pos] * valpiv;

            // Update trailing entries
            let jmax = if i <= nel2 { i } else { nel2 };
            for jj in 1..=jmax {
                let target = k1pos + jj;  // A[npiv+jj, col]
                A[target] -= A[k1pos] * A[apos + jj];
                // A[npiv+jj, npiv+i] -= L[npiv+i, npiv] * U[npiv, npiv+jj]
                // = (old_i/d) * old_j = rank-1 update
            }
        }

    } else {
        // === 2x2 PIVOT UPDATE ===
        let pospv1 = npiv + npiv * lda;         // A[npiv, npiv]
        let pospv2 = (npiv+1) + (npiv+1) * lda; // A[npiv+1, npiv+1]
        let offdag_old = npiv + (npiv+1) * lda;  // A[npiv, npiv+1] (upper)
        let offdag = (npiv+1) + npiv * lda;       // A[npiv+1, npiv] (lower)

        let d11 = A[pospv1];
        let d22 = A[pospv2];
        let detpiv = A[offdag];  // stored by FAC_I_LDLT
        let d12_orig = A[offdag_old];  // original off-diagonal

        // D^{-1} computation
        let a22 = d11 / detpiv;         // D^{-1}_{22}  (note: NOT d22/det)
        let a11 = d22 / detpiv;         // D^{-1}_{11}  (note: swapped!)
        let a12 = -d12_orig / detpiv;   // D^{-1}_{12}

        // Fix storage: move off-diagonal to lower triangle
        A[offdag] = d12_orig;       // A[npiv+1, npiv] = d12
        A[offdag_old] = 0.0;        // A[npiv, npiv+1] = 0

        // Process trailing rows (standard CPU path)
        // LPOS1 = A[npiv, npiv+2]  (first L row at first trailing column)
        // LPOS2 = A[npiv+1, npiv+2]

        // Within inner block (triangular): j2 = 1..nel2
        // jj starts at A[npiv, npiv+2] = pospv2 + lda - 1 in 1-based
        // In 0-based: jj_base = npiv + (npiv+2) * lda

        let mut jj = npiv + (npiv + 2) * lda;  // A[npiv, npiv+2]
        let mut ibeg = jj + 2;   // A[npiv+2, npiv+2] -- start of triangular update
        let mut iend = ibeg;     // initially only 1 entry

        for j2 in 0..nel2 {
            let col = npiv + 2 + j2;  // 0-based column
            let k1 = jj;              // A[npiv, col]
            let k2 = jj + 1;          // A[npiv+1, col]

            // Compute L*D^{-1} entries (negated)
            let mult1 = -(a11 * A[k1] + a12 * A[k2]);
            let mult2 = -(a12 * A[k1] + a22 * A[k2]);

            // Save originals to U storage (pivot rows)
            A[pospv1 + 2 + j2] = A[k1];  // U row 1
            A[pospv2 + 1 + j2] = A[k2];  // U row 2

            // Rank-2 Schur complement update (triangular)
            let mut uk1 = pospv1 + 2;     // start of saved U row 1
            let mut uk2 = pospv2 + 1;     // start of saved U row 2
            for irow in ibeg..=iend {
                A[irow] += mult1 * A[uk1] + mult2 * A[uk2];
                uk1 += 1;
                uk2 += 1;
            }

            // Store scaled L entries
            A[jj] = -mult1;       // L*D^{-1} entry for row npiv
            A[jj + 1] = -mult2;   // L*D^{-1} entry for row npiv+1

            jj += lda;
            ibeg += lda;
            iend += lda + 1;  // triangular region grows by 1 each iteration
        }

        // Beyond inner block (rectangular): rows iend_block..last_row-1
        iend -= 1;  // adjust after loop
        for j2 in 0..(last_row - iend_block) {
            let row_shift = j2 * lda;
            let jj_loc = jj + row_shift;
            let ibeg_loc = ibeg + row_shift;
            let iend_loc = iend + row_shift;

            let k1 = jj_loc;
            let k2 = jj_loc + 1;

            let mult1 = -(a11 * A[k1] + a12 * A[k2]);
            let mult2 = -(a12 * A[k1] + a22 * A[k2]);

            A[pospv1 + 2 + nel2 + j2] = A[k1];
            A[pospv2 + 1 + nel2 + j2] = A[k2];

            let mut uk1 = pospv1 + 2;
            let mut uk2 = pospv2 + 1;
            for irow in ibeg_loc..=iend_loc {
                A[irow] += mult1 * A[uk1] + mult2 * A[uk2];
                uk1 += 1;
                uk2 += 1;
            }

            A[jj_loc] = -mult1;
            A[jj_loc + 1] = -mult2;
        }
    }

    return ifinb;
}
```

### between_block_update_ldlt (corresponds to FAC_SQ_LDLT, GEMM-only call)

```rust
fn between_block_update_ldlt(
    A: &mut [f64], lda: usize,
    ibeg_block: usize,   // 0-based start of pivot block
    iend_block: usize,   // 0-based one-past-end of pivot block
    npiv: usize,         // current pivot count
    last_col_gemm: usize,  // update columns up to here
    last_row_gemm: usize,  // update rows up to here
) {
    let npiv_block = npiv - ibeg_block;
    if npiv_block == 0 { return; }
    let nel1 = last_col_gemm - iend_block;
    if nel1 == 0 { return; }

    // Symmetric part: update [iend_block..last_col_gemm) x [iend_block..last_col_gemm)
    let blsize = ...; // blocking parameter
    for irow in (iend_block..last_col_gemm).step_by(blsize) {
        let block = min(blsize, last_col_gemm - irow);
        let upos = ibeg_block * lda + irow;     // U: A[irow, ibeg_block]
        let lpos = irow * lda + ibeg_block;      // L: A[ibeg_block, irow]
        let apos = irow * lda + irow;            // target: A[irow, irow]

        // C -= A * B   where A = U region (block x npiv_block), B = L region
        dgemm('N', 'N', block, last_col_gemm - irow, npiv_block,
              -1.0, &A[upos], lda,
              &A[lpos], lda,
              1.0, &mut A[apos], lda);
    }

    // Rectangular part: update [iend_block..last_col_gemm) x [last_col_gemm..last_row_gemm)
    if last_row_gemm > last_col_gemm {
        let upos = ibeg_block * lda + iend_block;
        let lpos = last_col_gemm * lda + ibeg_block;
        let apos = last_col_gemm * lda + iend_block;

        dgemm('N', 'N', nel1, last_row_gemm - last_col_gemm, npiv_block,
              -1.0, &A[upos], lda,
              &A[lpos], lda,
              1.0, &mut A[apos], lda);
    }
}
```

### inter_panel_update_ldlt (FAC_SQ_LDLT with TRSM, for inter-panel)

```rust
fn inter_panel_update_ldlt(
    A: &mut [f64], lda: usize, nass: usize,
    ibeg_blr: usize, iend_blr: usize, npiv: usize,
    last_row: usize,
) {
    let npiv_block = npiv - ibeg_blr;
    if npiv_block == 0 { return; }

    // For PIVOT_OPTION <= 1, a TRSM is needed to solve for L entries
    // beyond the BLR panel. For PIVOT_OPTION=3, the MQ update already
    // handled all rows, so only GEMM is needed.

    // For PIVOT_OPTION=3:
    // CALL_TRSM = false (PIVOT_OPTION <= 1 evaluates to false)
    // Just do GEMM for columns iend_blr..nass-1 and rows nass..nfront-1

    // The GEMM updates:
    //   columns iend_blr..nass using pivots ibeg_blr..npiv
    //   rows iend_blr..nass AND nass..nfront

    // Symmetric GEMM for [iend_blr..nass) x [iend_blr..nass)
    // Rectangular GEMM for [iend_blr..nass) x [nass..last_row)

    // ... same structure as between_block_update_ldlt ...
}
```

---

## CRITICAL DETAILS FOR CORRECTNESS

### 1. The "U" storage pattern

After factoring column `k`, the original (unscaled) values from the L column are saved in the PIVOT ROW:
- For 1x1 pivot at position `k`: `A[k+1, k], A[k+2, k], ...` (the original values before scaling) are copied to `A[k, k+1], A[k, k+2], ...` (the "U" row). The L column then stores `original / d`.
- For 2x2 pivot at positions `k, k+1`: The originals from rows `k` and `k+1` of trailing columns are saved to `A[k, k+2..], A[k+1, k+2..]` (the two "U" rows). The L entries then store `D^{-1} * original`.

This is NOT simply "the upper triangle stores U and lower stores L". The U row and L column share the same values conceptually (by symmetry of the original matrix), but they store DIFFERENT things after factorization:
- **U row (above diagonal)**: original unscaled values = `D * L^T`
- **L column (below diagonal)**: scaled values = `L * D^{-1}` (which is just `L` for 1x1, but `L * D^{-1}` for the within-panel update interpretation)

Wait, let me be very precise. For a 1x1 pivot with value `d` at position `k`:
- `L[i, k]` = `original_a[i, k] / d` — stored in column k, rows i > k
- `U[k, j]` = `original_a[k, j]` = `original_a[j, k]` (by symmetry before factoring) = `d * L[j, k]` — stored in row k, columns j > k

The Schur complement update is: `A[i,j] -= L[i,k] * U[k,j] = L[i,k] * d * L[j,k]`

This is mathematically `A -= L * D * L^T` applied one column at a time.

### 2. The DETPIV storage for 2x2 pivots

When a 2x2 pivot is selected in FAC_I_LDLT, `DETPIV = d11*d22 - d12^2` is stored at position `A[NPIV+2, NPIV+1]` (the sub-diagonal of the 2x2 block, 1-based). In 0-based: `A[(npiv+1) + npiv * lda]`.

This DETPIV value is then READ by FAC_MQ_LDLT (not recomputed). The D^{-1} is computed as:
```
A22_of_Dinv = d11 / det    // note: d11 goes to (2,2) position
A11_of_Dinv = d22 / det    // note: d22 goes to (1,1) position
A12_of_Dinv = -d12 / det
```

The notation is confusing because MUMPS swaps the names: variable `A22` gets `A[POSPV1]/DETPIV` which is `d11/det`, and variable `A11` gets `SWOP/DETPIV = d22/det`. This is because in the inverse of `[[d11, d12], [d12, d22]]`, the (1,1) element is `d22/det` and the (2,2) element is `d11/det`.

### 3. The off-diagonal storage fix

After the 2x2 update, the off-diagonal is stored ONLY in the lower triangle:
- `A[npiv+1, npiv]` = original d12 (the off-diagonal of D)
- `A[npiv, npiv+1]` = 0 (zeroed)

This is important for the solve phase, which reads D from the lower triangle.

### 4. 2x2 pivot marker in IW (integer workspace)

After a 2x2 pivot, the column index of the SECOND pivot column is NEGATED in the IW array:
```
IW[col_index_position_of_npiv+1] = -IW[col_index_position_of_npiv+1]
```

This allows the solve phase and COPY2U_SCALEL to detect 2x2 pivots by checking for negative column indices.

### 5. Negative eigenvalue counting for 2x2 pivots

The eigenvalues of a 2x2 block `[[d11, d12], [d12, d22]]` have signs determined by the trace and determinant:
- If `det < 0`: one positive, one negative eigenvalue -> `NNEG += 1`
- If `det > 0` and `d22 < 0`: both negative -> `NNEG += 2`
- If `det > 0` and `d22 >= 0`: both positive -> `NNEG += 0`

### 6. MAXFROMM optimization

When `PIVOT_OPTION >= 3` and `IS_MAX_USEFUL` is true (i.e., `UU != 0`), the MQ update computes `MAXFROMM` = the maximum absolute value in the FIRST trailing column (column `npiv+pivsiz+1`) after the update. This is passed to the next call of FAC_I_LDLT as a shortcut: if the diagonal of the next candidate is large enough relative to MAXFROMM, the full column scan can be skipped. This is an optimization only and does not affect correctness.

### 7. The scan regions in the 1x1 pivot test

For PIVOT_OPTION=3, the pivot threshold test checks:
```
|diagonal| >= UU * max(AMAX, RMAX)
```
where:
- **AMAX** = max off-diagonal in column IPIV within rows `[NPIV+1, IEND_BLOCK]` — scanned by walking down the column within the inner block AND across to other columns (by symmetry)
- **RMAX** = max off-diagonal in column IPIV in rows `[IEND_BLOCK+1, NFRONT - KEEP(253) - NVSCHUR]` — the rows outside the inner block

For the 2x2 test, RMAX is refined to exclude the partner entry, and TMAX is the maximum in the partner's column excluding the off-diagonal.

### 8. INOPV = 2 handling

When FAC_I_LDLT returns `INOPV = 2`, it means no pivot was found within `[NPIV+1, IEND_BLOCK]`. The driver breaks out of the pivot loop but continues the inner block loop, which will immediately advance `IEND_BLOCK` and try again with more columns. This is how delayed pivots are implicitly handled — they stay in the frontal matrix at their current position and are retried when the inner block window expands.

When `INOPV = 2` is returned, the between-block GEMM (FAC_SQ_LDLT) is still called to update the new columns that will be tried next time.

### 9. Static pivoting

When static mode is active and no pivot is found (`INOPV = 1`), the driver sets `INOPV = -1` and retries. FAC_I_LDLT then takes the diagonal at `(NPIV+1, NPIV+1)` unconditionally, replacing it with `±CSEUIL` if it's too small (where `CSEUIL = max(DKEEP(1), SEUIL)`).

### 10. What happens to delayed pivots at the end

When the while-loop over NASS terminates (either normally or with `LASTPANEL = true`), the number of successfully eliminated pivots is `NPIV`. The remaining `NASS - NPIV` columns that could not be pivoted are "delayed" — they are added to the parent node's frontal matrix as additional fully-summed variables. The contribution block (rows/columns `NPIV+1..NFRONT`) is passed up to the parent. The `NE_STEPS` and `ND_STEPS` arrays track how many variables were actually eliminated vs total.