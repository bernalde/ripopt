# Makefile for the GAMS ripopt solver link
#
# Usage:
#   make                       — build libGamsRipopt.dylib (macOS) or .so (Linux)
#   make GAMS_PATH=/path/gams  — override GAMS installation path
#   make clean                 — remove build artifacts (the shared library)
#   make clean-bench-<size>    — remove benchmark outputs for one size
#                                so `make bench-<size>` actually re-runs
#                                (sizes: smoke, small, medium, large)
#   make clean-bench           — remove all benchmark outputs
#   make clean-all             — clean + clean-bench
#
# Prerequisites:
#   - GAMS installation (for C API headers and stubs)
#   - libripopt built: cargo build --release  (in repo root)

# --- Paths ---

# GAMS system directory (contains apifiles/)
GAMS_PATH ?= /Library/Frameworks/GAMS.framework/Versions/Current/Resources

# ripopt repo root
RIPOPT_ROOT ?= ..

# GAMS C API source directory
GAMS_API = $(GAMS_PATH)/apifiles/C/api

# --- Platform detection ---

UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
  SHLIB_EXT   = dylib
  SHLIB_FLAGS = -dynamiclib -install_name @rpath/libGamsRipopt.dylib
  RPATH_FLAG  = -Wl,-rpath,$(RIPOPT_ROOT)/target/release
else
  SHLIB_EXT   = so
  SHLIB_FLAGS = -shared
  RPATH_FLAG  = -Wl,-rpath,'$$ORIGIN'
endif

# --- Compiler ---

CC       = cc
CFLAGS   = -O2 -Wall -Wextra -Wno-unused-parameter \
           -I$(RIPOPT_ROOT) -I$(GAMS_API)
LDFLAGS  = -L$(RIPOPT_ROOT)/target/release -lripopt $(RPATH_FLAG)

# --- Targets ---

TARGET = libGamsRipopt.$(SHLIB_EXT)

GAMS_SRCS = $(GAMS_API)/gmomcc.c $(GAMS_API)/gevmcc.c

.PHONY: all clean install test help

help:
	@echo "Targets:"
	@echo "  all                  Build libGamsRipopt.dylib/.so (default)"
	@echo "  install              Build and install library + register RIPOPT solver in GAMS"
	@echo "  test                 install + run test_hs071.gms"
	@echo "  clean                Remove build artifacts (shared library)"
	@echo "  clean-bench-smoke    Remove smoke benchmark outputs"
	@echo "  clean-bench-small    Remove small benchmark outputs"
	@echo "  clean-bench-medium   Remove medium benchmark outputs"
	@echo "  clean-bench-large    Remove large benchmark outputs"
	@echo "  clean-bench          Remove all benchmark outputs"
	@echo "  clean-all            clean + clean-bench"
	@echo "  bench-smoke          Run ~10 problems, TIMELIMIT=30s"
	@echo "  bench-small          Run 50 problems, TIMELIMIT=60s"
	@echo "  bench-medium         Run 50 problems, TIMELIMIT=300s"
	@echo "  bench-large          Run 50 problems, TIMELIMIT=900s"
	@echo "  bench-all            Run all four benchmark sizes"
	@echo "  testsets             Regenerate nlpbench/testsets/*.gms"
	@echo ""
	@echo "Variables:"
	@echo "  GAMS_PATH=<path>     Override GAMS installation directory"
	@echo "  TIMELIMIT=<sec>      Override per-instance timeout for bench targets"
	@echo "  VERSION=vX.Y.Z       Override version tag in report filenames"
	@echo "  REPORT_SOLVERS=...   Solvers to include in reports (default: ripopt ipopt)"

all: $(TARGET)

$(TARGET): gams_ripopt.c $(GAMS_SRCS)
	$(CC) $(CFLAGS) $(SHLIB_FLAGS) -o $@ $^ $(LDFLAGS)

install: $(TARGET)
	cp $(TARGET) $(GAMS_PATH)/
	cp $(RIPOPT_ROOT)/target/release/libripopt.$(SHLIB_EXT) $(GAMS_PATH)/
ifeq ($(UNAME),Darwin)
	install_name_tool -change \
		"$$(otool -L $(GAMS_PATH)/$(TARGET) | grep libripopt | awk '{print $$1}')" \
		@loader_path/libripopt.$(SHLIB_EXT) \
		$(GAMS_PATH)/$(TARGET)
	install_name_tool -id @loader_path/$(TARGET) $(GAMS_PATH)/$(TARGET)
endif
	@grep -q '^RIPOPT ' $(GAMS_PATH)/gmscmpun.txt 2>/dev/null || \
		(cp $(GAMS_PATH)/gmscmpun.txt $(GAMS_PATH)/gmscmpun.txt.bak && \
		awk '/^DEFAULTS$$/{print "RIPOPT 11 5 00010203040506070809 1 0 2 NLP DNLP RMINLP"; print "gmsgenus.run"; print "gmsgenux.out"; print "libGamsRipopt.$(SHLIB_EXT) rip 1 0"}1' \
		$(GAMS_PATH)/gmscmpun.txt.bak > $(GAMS_PATH)/gmscmpun.txt)
	@echo "Installed to $(GAMS_PATH)"

test: install
	$(GAMS_PATH)/gams test_hs071.gms

clean:
	rm -f $(TARGET)

# --- nlpbench benchmark targets ---
#
# Each bench-<size> target:
#   1. runs ripopt and ipopt on gams/testsets/<size>.gms via run_nlpbench.sh
#   2. generates BENCHMARK_REPORT_<size>_<version>.md in gams/
#
# Per-instance timeout is tuned per size bucket. Override with TIMELIMIT=...
# Override version with VERSION=vX.Y.Z (default reads Cargo.toml).

REPORT_SOLVERS ?= ripopt ipopt
VERSION ?= $(shell grep '^version' $(RIPOPT_ROOT)/Cargo.toml | head -1 | sed 's/.*"\(.*\)".*/v\1/')

.PHONY: bench-smoke bench-small bench-medium bench-large bench-all testsets

bench-smoke:  TIMELIMIT ?= 30
bench-small:  TIMELIMIT ?= 60
bench-medium: TIMELIMIT ?= 300
bench-large:  TIMELIMIT ?= 900

bench-smoke bench-small bench-medium bench-large:
	@set -e; \
	size=$(patsubst bench-%,%,$@); \
	echo "=== bench $$size (TIMELIMIT=$(TIMELIMIT), version=$(VERSION)) ==="; \
	TIMELIMIT=$(TIMELIMIT) bash nlpbench/run_nlpbench.sh $$size.gms $(REPORT_SOLVERS); \
	traces=""; \
	for s in $(REPORT_SOLVERS); do \
	    trace=nlpbench/runsolver/$$size.gms_$$s.csv; \
	    if [ ! -f "$$trace" ]; then \
	        echo "error: expected trace $$trace not produced" >&2; \
	        exit 1; \
	    fi; \
	    traces="$$traces --trace $$trace"; \
	done; \
	python3 nlpbench/nlpbench_report.py $$traces \
	    --out nlpbench/BENCHMARK_REPORT_$${size}_$(VERSION).md \
	    --version $(VERSION); \
	echo "report: gams/nlpbench/BENCHMARK_REPORT_$${size}_$(VERSION).md"

bench-all: bench-smoke bench-small bench-medium bench-large

# Regenerate nlpbench/testsets/*.gms from the upstream conopt3 reference trace.
testsets:
	python3 nlpbench/gen_testsets.py

# --- Benchmark cleanup ---
#
# The upstream runsolver/Makefile caches per-instance trace CSVs at
#   nlpbench/runsolver/<solver>/<instance>.csv
# and the concatenated trace at
#   nlpbench/runsolver/<testset>.gms_<solver>.csv
# Make's dependency tracking means a rerun skips instances whose CSV already
# exists. To actually re-run, those files must be removed first.
#
# clean-bench-<size> reads nlpbench/testsets/<size>.gms and removes:
#   - the concatenated trace for each solver in $(REPORT_SOLVERS)
#   - the per-instance .csv/.log/.lst for each instance and each solver
#   - the BENCHMARK_REPORT_<size>_*.md reports
# Per-instance files are shared across testsets, so cleaning `small` does NOT
# clean instances that only `medium` or `large` cover.

.PHONY: clean-bench clean-bench-smoke clean-bench-small clean-bench-medium clean-bench-large clean-all

clean-bench-smoke clean-bench-small clean-bench-medium clean-bench-large:
	@size=$(patsubst clean-bench-%,%,$@); \
	testset="nlpbench/testsets/$${size}.gms"; \
	echo "=== cleaning $$size benchmark outputs ==="; \
	for s in $(REPORT_SOLVERS); do \
	    rm -fv "nlpbench/runsolver/$${size}.gms_$${s}.csv" 2>/dev/null || true; \
	    if [ -f "$$testset" ]; then \
	        while IFS= read -r inst; do \
	            case "$$inst" in ''|\#*) continue ;; esac; \
	            rm -f "nlpbench/runsolver/$${s}/$${inst}.csv" \
	                  "nlpbench/runsolver/$${s}/$${inst}.log" \
	                  "nlpbench/runsolver/$${s}/$${inst}.lst"; \
	        done < "$$testset"; \
	    else \
	        echo "warning: $$testset not found; skipping per-instance cleanup" >&2; \
	    fi; \
	done; \
	rm -fv nlpbench/BENCHMARK_REPORT_$${size}_*.md 2>/dev/null || true

clean-bench: clean-bench-smoke clean-bench-small clean-bench-medium clean-bench-large

clean-all: clean clean-bench
