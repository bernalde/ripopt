#!/bin/bash
# install.sh — Register ripopt as a GAMS solver.
#
# Usage:
#   ./install.sh                       # uses default GAMS path
#   GAMS_PATH=/path/to/gams ./install.sh
#
# What it does:
#   1. Copies libGamsRipopt and libripopt into the GAMS system directory
#   2. Adds a RIPOPT entry to gmscmpun.txt (the solver registration file)

set -euo pipefail

GAMS_PATH="${GAMS_PATH:-/Library/Frameworks/GAMS.framework/Versions/Current/Resources}"
RIPOPT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"

UNAME="$(uname -s)"
if [ "$UNAME" = "Darwin" ]; then
    EXT="dylib"
else
    EXT="so"
fi

LINK_LIB="$RIPOPT_ROOT/gams/libGamsRipopt.$EXT"
SOLVER_LIB="$RIPOPT_ROOT/target/release/libripopt.$EXT"
CMPTXT="$GAMS_PATH/gmscmpun.txt"

# --- Checks ---

if [ ! -f "$LINK_LIB" ]; then
    echo "Error: $LINK_LIB not found. Run 'make -C gams' first." >&2
    exit 1
fi

if [ ! -f "$SOLVER_LIB" ]; then
    echo "Error: $SOLVER_LIB not found. Run 'cargo build --release' first." >&2
    exit 1
fi

if [ ! -d "$GAMS_PATH" ]; then
    echo "Error: GAMS directory $GAMS_PATH not found." >&2
    exit 1
fi

# --- Copy libraries ---

echo "Copying libGamsRipopt.$EXT → $GAMS_PATH/"
cp "$LINK_LIB" "$GAMS_PATH/"

echo "Copying libripopt.$EXT → $GAMS_PATH/"
cp "$SOLVER_LIB" "$GAMS_PATH/"

# Fix library paths on macOS so libGamsRipopt finds libripopt in the same dir
if [ "$UNAME" = "Darwin" ]; then
    # Find the actual recorded path for libripopt and rewrite it
    OLD_PATH=$(otool -L "$GAMS_PATH/libGamsRipopt.dylib" | grep libripopt | awk '{print $1}')
    if [ -n "$OLD_PATH" ]; then
        echo "Rewriting libripopt path: $OLD_PATH → @loader_path/libripopt.dylib"
        install_name_tool -change \
            "$OLD_PATH" \
            "@loader_path/libripopt.dylib" \
            "$GAMS_PATH/libGamsRipopt.dylib"
    fi

    install_name_tool -id \
        "@loader_path/libGamsRipopt.dylib" \
        "$GAMS_PATH/libGamsRipopt.dylib" 2>/dev/null || true
fi

# --- Register solver in gmscmpun.txt ---

if grep -q "^RIPOPT " "$CMPTXT" 2>/dev/null; then
    echo "RIPOPT already registered in $CMPTXT"
else
    echo "Adding RIPOPT to $CMPTXT (before DEFAULTS section)"
    # Insert RIPOPT entry before the DEFAULTS line
    sed -i.bak '/^DEFAULTS$/i\
RIPOPT 11 5 00010203040506070809 1 0 2 NLP DNLP RMINLP\
gmsgenus.run\
gmsgenux.out\
libGamsRipopt.dylib rip 1 0\
' "$CMPTXT"
    echo "Done."
fi

echo ""
echo "Installation complete. Test with:"
echo "  gams gams/test_hs071.gms"
