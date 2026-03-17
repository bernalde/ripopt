SHELL = /bin/bash
.PHONY: help test test-c test-iterative install uninstall clean cutest-install cutest-prepare cutest-run cutest-full cutest-report cutest cutest-maxiter cutest-smoke cutest-large hs-run large-scale electrolyte-run opf-run cho-run benchmark benchmark-report

# Show available targets
help:
	@echo "Usage: make <target> [CUTEST_MAX_N=N]"
	@echo ""
	@echo "Build & Install:"
	@echo "  install          Build and install ripopt binary, shared library, and Pyomo plugin"
	@echo "  uninstall        Remove ripopt binary, shared library, and Pyomo plugin"
	@echo "  clean            Remove build artifacts (cargo clean)"
	@echo ""
	@echo "Testing:"
	@echo "  test             Run unit/integration tests (cargo test)"
	@echo "  test-iterative   Run iterative/hybrid solver tests (release, large-scale)"
	@echo ""
	@echo "Benchmarks:"
	@echo "  benchmark        Full benchmark: HS + CUTEst + domain + large-scale + report"
	@echo "  benchmark-report Generate unified report from existing results"
	@echo "  hs-run           Run HS suite benchmark (ripopt + ipopt)"
	@echo "  electrolyte-run  Run electrolyte thermodynamics benchmark"
	@echo "  opf-run          Run AC optimal power flow benchmark"
	@echo "  cho-run          Run CHO parameter estimation benchmark"
	@echo "  large-scale      Run synthetic large-scale tests (up to 100K vars)"
	@echo "  test-c           Build shared library and run C API examples"
	@echo ""
	@echo "CUTEst Benchmarks:"
	@echo "  cutest-install   Install CUTEst framework to ~/.local/cutest"
	@echo "  cutest           Full pipeline: prepare, run, report"
	@echo "  cutest-prepare   Compile SIF problems -> .dylib (reads problem_list.txt)"
	@echo "  cutest-run       Run benchmark on all prepared problems (results.json)"
	@echo "  cutest-report    Generate comparison report from results.json"
	@echo "  cutest-full      Full benchmark: all 1542 SIF problems (full_results.json)"
	@echo "  cutest-smoke     Quick smoke test (~20 small problems, n+m < 50)"
	@echo "  cutest-large     Large problems (n+m >= 100, sparse solver)"
	@echo "  cutest-maxiter   Re-run the 35 MaxIterations failures"
	@echo ""
	@echo "Environment variables:"
	@echo "  CUTEST_MAX_N=N   Skip problems with n > N (default: 100)"

# Detect shared library extension
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  DYLIB_EXT := dylib
else
  DYLIB_EXT := so
endif

# Install ripopt binary and shared library
install:
	cargo build --release
	cargo install --path . --bin ripopt
	mkdir -p ~/.local/lib
	cp target/release/libripopt.$(DYLIB_EXT) ~/.local/lib/
	@echo ""
	@echo "Installed:"
	@echo "  ripopt binary      -> ~/.cargo/bin/ripopt"
	@echo "  libripopt.$(DYLIB_EXT)  -> ~/.local/lib/libripopt.$(DYLIB_EXT)"
	@echo ""
	@if echo "$$PATH" | tr ':' '\n' | grep -qx "$$HOME/.cargo/bin"; then \
		echo "Verify: ripopt --version"; \
	else \
		echo "WARNING: ~/.cargo/bin is not on your PATH."; \
		echo "Add it by appending this to your shell profile (~/.bashrc, ~/.zshrc, etc.):"; \
		echo "  export PATH=\"\$$HOME/.cargo/bin:\$$PATH\""; \
		echo "Then restart your shell or run: source ~/.bashrc"; \
	fi
	@if command -v pip >/dev/null 2>&1; then \
		echo ""; \
		echo "Installing Pyomo solver plugin..."; \
		pip install ./pyomo-ripopt && \
		echo "  pyomo-ripopt       -> installed via pip"; \
	else \
		echo ""; \
		echo "NOTE: pip not found. To install the Pyomo solver plugin later, run:"; \
		echo "  pip install ./pyomo-ripopt"; \
	fi
	@if echo "$$LD_LIBRARY_PATH:$$DYLD_LIBRARY_PATH" | tr ':' '\n' | grep -qx "$$HOME/.local/lib"; then \
		true; \
	else \
		echo ""; \
		echo "NOTE: To use the shared library, ensure ~/.local/lib is in your library path:"; \
		echo "  export LD_LIBRARY_PATH=\"\$$HOME/.local/lib:\$$LD_LIBRARY_PATH\""; \
	fi

# Uninstall ripopt binary and shared library
uninstall:
	cargo uninstall ripopt 2>/dev/null || true
	pip uninstall -y pyomo-ripopt 2>/dev/null || true
	rm -f ~/.local/lib/libripopt.$(DYLIB_EXT)
	@echo "Uninstalled ripopt binary and shared library."

# Remove build artifacts
clean:
	cargo clean

# Run unit/integration tests (including C API tests via Rust FFI)
test:
	cargo test

# Run iterative/hybrid linear solver integration tests (includes ignored large-scale tests)
test-iterative:
	cargo test --release -- --ignored iterative hybrid

# Build shared library and compile/run all C API examples
test-c:
	cargo build --release
	@echo "=== Compiling and running C API examples ==="
	cc examples/c_api_test.c -I. -Ltarget/release -lripopt \
		-Wl,-rpath,$(CURDIR)/target/release -o target/release/c_api_test -lm
	target/release/c_api_test
	cc examples/c_rosenbrock.c -I. -Ltarget/release -lripopt \
		-Wl,-rpath,$(CURDIR)/target/release -o target/release/c_rosenbrock -lm
	target/release/c_rosenbrock
	cc examples/c_hs035.c -I. -Ltarget/release -lripopt \
		-Wl,-rpath,$(CURDIR)/target/release -o target/release/c_hs035 -lm
	target/release/c_hs035
	cc examples/c_example_with_options.c -I. -Ltarget/release -lripopt \
		-Wl,-rpath,$(CURDIR)/target/release -o target/release/c_example_with_options -lm
	target/release/c_example_with_options
	@echo "=== All C API examples passed ==="

# Install CUTEst framework to ~/.local/cutest
cutest-install:
	bash install_cutest.sh

# Prepare CUTEst problems (compile SIF -> .dylib)
cutest-prepare:
	source ~/.local/cutest/env.sh && bash cutest_suite/prepare.sh

# Run CUTEst benchmark (saves results to cutest_suite/results.json)
cutest-run:
	RESULTS_FILE=cutest_suite/results.json \
	cargo run --bin cutest_suite --features cutest,ipopt-native --release \
		2> >(tee cutest_suite/benchmark_stderr.txt >&2)

# Run full CUTEst benchmark (all 1542 SIF problems)
cutest-full:
	PROBLEM_LIST=cutest_suite/problem_list_full.txt \
	RESULTS_FILE=cutest_suite/full_results.json \
	cargo run --bin cutest_suite --features cutest,ipopt-native --release \
		2> >(tee cutest_suite/full_stderr.txt >&2)

# Generate comparison report from results
cutest-report:
	python cutest_suite/compare.py cutest_suite/results.json

# Smoke test: ~20 fast, diverse problems (unconstrained, bounds, equality,
# inequality, NE-to-LS, mixed). All n+m < 50, finishes in seconds.
cutest-smoke:
	RESULTS_FILE=cutest_suite/smoke_results.json \
	cargo run --bin cutest_suite --features cutest,ipopt-native --release -- \
		ROSENBR BEALE CUBE DENSCHNB BROWNBS \
		HS6 HS10 HS35 HS71 HS106 \
		MARATOS BT1 HS13 HS57 S316-322 \
		BEALENE BROWNBSNE ENGVAL2NE \
		PALMER1D PALMER5D \
		2> >(tee cutest_suite/smoke_stderr.txt >&2)
	@echo "Smoke test complete."

# Large problems: n+m >= 100, exercises the sparse multifrontal solver.
# Includes over-constrained (m >> n), large-n, and NE systems.
cutest-large:
	RESULTS_FILE=cutest_suite/large_results.json \
	cargo run --bin cutest_suite --features cutest,ipopt-native --release -- \
		HIMMELBI AIRPORT LAKES SWOPF QPCBLEND QPNBLEND SPANHYD LINSPANH \
		ACOPP14 ACOPR14 ACOPP30 ACOPR30 BATCH CORE1 MSS1 HAIFAM FEEDLOC \
		DISCS DECONVBNE DECONVNE NET1 HYDCAR20 \
		TFI1 TFI2 TFI3 EXPFITB ELATTAR CRESC50 EXPFITC \
		DUALC1 DUALC2 DUALC5 DUALC8 \
		PT HET-Z OET1 OET2 OET3 OET4 KSIP \
		SIPOW1 SIPOW1M SIPOW2 SIPOW2M SIPOW3 SIPOW4 \
		TAXR13322 GOFFIN \
		2> >(tee cutest_suite/large_stderr.txt >&2)
	@echo "Large problem benchmark complete."

# Run only the 35 MaxIterations failures for fast iteration
cutest-maxiter:
	RESULTS_FILE=cutest_suite/maxiter_results.json \
	cargo run --bin cutest_suite --features cutest,ipopt-native --release -- \
		AIRCRFTB ALLINIT BENNETT5LS BIGGS3 BIGGS5 BOX2 BQPGABIM CHWIRUT2LS \
		DJTL HAHN1LS MEYER3 MGH10LS MINSURF OSBORNEA RAT43LS SIM2BQP THURBERLS \
		AIRCRFTA CONCON DECONVC DECONVNE DECONVU DNIEPER HS109 HS67 HS83 HS84 \
		HS85 HS99EXP MCONCON OPTCNTRL QC QCNEW SSINE TRIGGER \
		2> >(tee cutest_suite/maxiter_stderr.txt >&2)

# Full CUTEst pipeline: prepare, run, report
cutest: cutest-prepare cutest-run cutest-report
	@echo "CUTEst benchmark complete. Report: cutest_suite/CUTEST_REPORT.md"

# Run HS suite benchmark (ripopt + ipopt native)
hs-run:
	cargo run --bin hs_suite --features hs,ipopt-native --release \
		> hs_suite/hs_ripopt_results.json \
		2> >(tee hs_suite/hs_ripopt_stderr.txt >&2)
	cargo run --bin ipopt_native --features ipopt-native --release \
		> hs_suite/hs_ipopt_native_results.json \
		2> >(tee hs_suite/hs_ipopt_native_stderr.txt >&2)
	@echo "HS benchmark complete."

# Run synthetic large-scale tests (5K to 100K variables)
large-scale:
	cargo test --release --features ipopt-native -- --ignored large_scale_vs_ipopt --nocapture \
		2>&1 | tee large_scale_results.txt
	@echo "Large-scale benchmark complete. Output: large_scale_results.txt"

# Run electrolyte thermodynamics benchmark (ripopt + ipopt)
electrolyte-run:
	RESULTS_FILE=electrolyte_results.json \
	cargo run --release --features ipopt-native --example electrolyte_benchmark \
		2> >(tee electrolyte_stderr.txt >&2)
	@echo "Electrolyte benchmark complete."

# Run AC optimal power flow benchmark (ripopt + ipopt)
opf-run:
	RESULTS_FILE=opf_results.json \
	cargo run --release --features ipopt-native --example opf_benchmark \
		2> >(tee opf_stderr.txt >&2)
	@echo "OPF benchmark complete."

# Run CHO parameter estimation benchmark (ripopt + ipopt)
cho-run:
	RESULTS_FILE=cho_results.json \
	cargo run --release --features ipopt-native --example cho_benchmark \
		2> >(tee cho_stderr.txt >&2)
	@echo "CHO benchmark complete."

# Full benchmark: HS + CUTEst + domain + large-scale + unified report
benchmark: hs-run cutest-run electrolyte-run opf-run cho-run large-scale benchmark-report
	@echo "Full benchmark complete. Report: BENCHMARK_REPORT.md"

# Generate unified benchmark report from existing results
benchmark-report:
	python benchmark_report.py
	@echo "Report: BENCHMARK_REPORT.md"
