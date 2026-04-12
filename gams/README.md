# GAMS solver link for ripopt

This directory builds `libGamsRipopt`, a shared library that registers ripopt
as a GAMS NLP solver. Once installed, a GAMS model can invoke ripopt with

```
option nlp = ripopt;
solve mymodel using nlp minimizing obj;
```

## Files

- `gams_ripopt.c` — the solver link. Translates between the GAMS Modeling
  Object (GMO) API and ripopt's C API (`ripopt.h`). Entry points are
  `ripCreate`, `ripFree`, `ripReadyAPI`, and `ripCallSolver`.
- `Makefile` — builds `libGamsRipopt.{dylib,so}` and installs it into a GAMS
  installation.
- `install.sh` — convenience wrapper around `make install` that auto-detects
  the GAMS path on macOS and Linux.
- `test_hs071.gms` — GAMS model of Hock-Schittkowski problem 71. Used by
  `make test` to verify the solver link end-to-end.

## Prerequisites

- A working GAMS installation (the Makefile auto-detects
  `/Library/Frameworks/GAMS.framework/...` on macOS; override with
  `GAMS_PATH=/path/to/gams`).
- `libripopt` built in the repo root:

  ```
  cargo build --release
  ```

## Build

```
make -C gams
```

Produces `gams/libGamsRipopt.{dylib,so}`, linked against
`../target/release/libripopt.{dylib,so}`.

## Install

```
sudo make -C gams install
```

Copies `libGamsRipopt` and `libripopt` into the GAMS system directory and
registers a `RIPOPT` entry in `gmscmpun.txt` so GAMS sees ripopt as an
available NLP solver. On macOS, rewrites the install-name of `libripopt`
inside `libGamsRipopt` to `@loader_path/libripopt.dylib` so the loader
resolves it from the GAMS directory.

## Test

```
sudo make -C gams test
```

Runs `test_hs071.gms` through GAMS. The model aborts on an objective
mismatch (`|obj - 17.014| > 1e-2`), an unexpected solve status, or an
unexpected model status — so a clean run is a strong end-to-end check.

## Capabilities

Registered model types: `NLP`, `DNLP`, `RMINLP`. Mixed-integer and conic
problem types are not supported by the underlying ripopt solver.
