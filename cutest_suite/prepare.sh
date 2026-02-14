#!/usr/bin/env bash
# Compile CUTEst SIF problems into shared libraries (.dylib) for the benchmark suite.
#
# Usage:
#   bash prepare.sh ROSENBR HS35 HS71        # specific problems
#   bash prepare.sh --from-file problem_list.txt  # from file
#   bash prepare.sh                           # uses default problem_list.txt
#
# Requires: source ~/.local/cutest/env.sh (sets CUTEST, SIFDECODE, MASTSIF)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROBLEMS_DIR="${SCRIPT_DIR}/problems"

# Source CUTEst environment
if [[ -f ~/.local/cutest/env.sh ]]; then
    source ~/.local/cutest/env.sh
else
    echo "ERROR: ~/.local/cutest/env.sh not found. Install CUTEst first."
    exit 1
fi

# Verify required tools
for cmd in sifdecoder gfortran; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "ERROR: $cmd not found in PATH"
        exit 1
    fi
done

if [[ -z "${MASTSIF:-}" ]]; then
    echo "ERROR: MASTSIF not set. Source ~/.local/cutest/env.sh"
    exit 1
fi

# Parse arguments
PROBLEM_NAMES=()
if [[ $# -eq 0 ]]; then
    # Default: read from problem_list.txt
    LIST_FILE="${SCRIPT_DIR}/problem_list.txt"
    if [[ ! -f "$LIST_FILE" ]]; then
        echo "ERROR: No arguments and problem_list.txt not found"
        exit 1
    fi
    while IFS= read -r line; do
        line="${line%%#*}"  # strip comments
        line="$(echo "$line" | xargs)"  # trim whitespace
        [[ -n "$line" ]] && PROBLEM_NAMES+=("$line")
    done < "$LIST_FILE"
elif [[ "$1" == "--from-file" ]]; then
    LIST_FILE="$2"
    while IFS= read -r line; do
        line="${line%%#*}"
        line="$(echo "$line" | xargs)"
        [[ -n "$line" ]] && PROBLEM_NAMES+=("$line")
    done < "$LIST_FILE"
else
    PROBLEM_NAMES=("$@")
fi

mkdir -p "$PROBLEMS_DIR"

echo "Preparing ${#PROBLEM_NAMES[@]} CUTEst problems..."

SUCCESS=0
FAIL=0
SKIP=0

for NAME in "${PROBLEM_NAMES[@]}"; do
    SIF_FILE="${MASTSIF}/${NAME}.SIF"
    if [[ ! -f "$SIF_FILE" ]]; then
        echo "  SKIP $NAME (${NAME}.SIF not found in MASTSIF)"
        ((SKIP++)) || true
        continue
    fi

    # Check if already prepared
    if [[ -f "${PROBLEMS_DIR}/lib${NAME}.dylib" && -f "${PROBLEMS_DIR}/${NAME}_OUTSDIF.d" ]]; then
        echo "  OK   $NAME (already prepared)"
        ((SUCCESS++)) || true
        continue
    fi

    # Create temporary working directory
    WORK_TMPDIR="$(mktemp -d)"

    echo -n "  BUILD $NAME ... "

    # Run sifdecoder
    if ! (cd "$WORK_TMPDIR" && sifdecoder "$SIF_FILE" > /dev/null 2>&1); then
        echo "FAILED (sifdecoder)"
        rm -rf "$WORK_TMPDIR"
        ((FAIL++)) || true
        continue
    fi

    # Check that required files were generated
    if [[ ! -f "$WORK_TMPDIR/ELFUN.f" || ! -f "$WORK_TMPDIR/GROUP.f" || ! -f "$WORK_TMPDIR/RANGE.f" ]]; then
        echo "FAILED (missing Fortran files)"
        rm -rf "$WORK_TMPDIR"
        ((FAIL++)) || true
        continue
    fi

    # Collect all .f files (ELFUN.f, GROUP.f, RANGE.f, possibly EXTER.f, SETTYP.f)
    FORTRAN_FILES=("$WORK_TMPDIR"/*.f)

    # Compile into shared library
    if ! gfortran -shared -fPIC -O2 -o "$WORK_TMPDIR/lib${NAME}.dylib" "${FORTRAN_FILES[@]}" 2>/dev/null; then
        echo "FAILED (gfortran)"
        rm -rf "$WORK_TMPDIR"
        ((FAIL++)) || true
        continue
    fi

    # Copy outputs
    cp "$WORK_TMPDIR/lib${NAME}.dylib" "$PROBLEMS_DIR/"
    cp "$WORK_TMPDIR/OUTSDIF.d" "$PROBLEMS_DIR/${NAME}_OUTSDIF.d"

    # Cleanup
    rm -rf "$WORK_TMPDIR"

    echo "OK"
    ((SUCCESS++)) || true
done

echo ""
echo "Done: $SUCCESS prepared, $FAIL failed, $SKIP skipped"
