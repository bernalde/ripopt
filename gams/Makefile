# Makefile for the GAMS ripopt solver link
#
# Usage:
#   make                       — build libGamsRipopt.dylib (macOS) or .so (Linux)
#   make GAMS_PATH=/path/gams  — override GAMS installation path
#   make clean                 — remove build artifacts
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

.PHONY: all clean install test

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
