#!/usr/bin/env python3
"""
Unified benchmark report for ripopt vs Ipopt.

Reads results from:
  - hs_suite/hs_ripopt_results.json      (HS suite, ripopt)
  - hs_suite/hs_ipopt_native_results.json (HS suite, ipopt)
  - cutest_suite/results.json             (CUTEst 727 suite)

Produces a single BENCHMARK_REPORT.md with per-suite and combined statistics.

Usage:
    python benchmark_report.py [--output BENCHMARK_REPORT.md]
    python benchmark_report.py --baseline old_report.json  # regression detection
"""

import json
import math
import os
import sys
from collections import defaultdict
from datetime import datetime

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


# ---- Helpers ----

def is_solved(status):
    return status in ('Optimal', 'Acceptable')


def obj_diff(ro, co):
    """Relative objective difference with floor of 1.0."""
    if ro is None or co is None:
        return float('nan')
    if not isinstance(ro, (int, float)) or not isinstance(co, (int, float)):
        return float('nan')
    if math.isnan(ro) or math.isnan(co):
        return float('nan')
    denom = max(abs(co), abs(ro), 1.0)
    return abs(ro - co) / denom


def fmt_time(t):
    if t is None or (isinstance(t, float) and math.isnan(t)):
        return "N/A"
    if t >= 1.0:
        return f"{t:.2f}s"
    elif t >= 0.001:
        return f"{t*1000:.1f}ms"
    else:
        return f"{t*1e6:.0f}us"


def geo_mean(values):
    """Geometric mean of positive values."""
    pos = [v for v in values if v > 0]
    if not pos:
        return float('nan')
    return math.exp(sum(math.log(v) for v in pos) / len(pos))


def median(values):
    if not values:
        return float('nan')
    s = sorted(values)
    return s[len(s) // 2]


def compute_stats(diffs):
    if not diffs:
        return float('nan'), float('nan'), float('nan')
    return sum(diffs) / len(diffs), median(diffs), max(diffs)


# ---- Load results ----

def load_hs_results():
    """Load HS suite results (separate files for each solver)."""
    ripopt_path = os.path.join(SCRIPT_DIR, 'hs_suite', 'hs_ripopt_results.json')
    ipopt_path = os.path.join(SCRIPT_DIR, 'hs_suite', 'hs_ipopt_native_results.json')

    if not os.path.exists(ripopt_path) or not os.path.exists(ipopt_path):
        return None

    with open(ripopt_path) as f:
        ripopt_data = json.load(f)
    with open(ipopt_path) as f:
        ipopt_data = json.load(f)

    ripopt_by_num = {r['number']: r for r in ripopt_data}
    ipopt_by_num = {r['number']: r for r in ipopt_data}

    comparisons = []
    for num in sorted(set(ripopt_by_num.keys()) | set(ipopt_by_num.keys())):
        rr = ripopt_by_num.get(num, {})
        cr = ipopt_by_num.get(num, {})

        r_solved = is_solved(rr.get('status', ''))
        c_solved = is_solved(cr.get('status', ''))
        both = r_solved and c_solved
        od = obj_diff(rr.get('objective'), cr.get('objective')) if both else float('nan')

        comparisons.append({
            'name': f"HS{num:03d}",
            'suite': 'HS',
            'n': rr.get('n', cr.get('n', 0)),
            'm': rr.get('m', cr.get('m', 0)),
            'ripopt_status': rr.get('status', 'N/A'),
            'ipopt_status': cr.get('status', 'N/A'),
            'ripopt_obj': rr.get('objective', float('nan')),
            'ipopt_obj': cr.get('objective', float('nan')),
            'obj_diff': od,
            'ripopt_iters': rr.get('iterations', 0),
            'ipopt_iters': cr.get('iterations', 0),
            'ripopt_time': rr.get('solve_time', 0),
            'ipopt_time': cr.get('solve_time', 0),
            'ripopt_solved': r_solved,
            'ipopt_solved': c_solved,
            'both_solved': both,
            'passed': both and not math.isnan(od) and od < 1e-4,
        })

    return comparisons


def load_cutest_results(path=None):
    """Load CUTEst results (single file with solver field)."""
    if path is None:
        path = os.path.join(SCRIPT_DIR, 'cutest_suite', 'results.json')

    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None

    with open(path) as f:
        data = json.load(f)

    ripopt_by_name = {}
    ipopt_by_name = {}
    for r in data:
        if r['solver'] == 'ripopt':
            ripopt_by_name[r['name']] = r
        elif r['solver'] == 'ipopt':
            ipopt_by_name[r['name']] = r

    comparisons = []
    for name in sorted(set(ripopt_by_name.keys()) | set(ipopt_by_name.keys())):
        rr = ripopt_by_name.get(name, {})
        cr = ipopt_by_name.get(name, {})

        r_solved = is_solved(rr.get('status', ''))
        c_solved = is_solved(cr.get('status', ''))
        both = r_solved and c_solved
        od = obj_diff(rr.get('objective'), cr.get('objective')) if both else float('nan')

        comparisons.append({
            'name': name,
            'suite': 'CUTEst',
            'n': rr.get('n', cr.get('n', 0)),
            'm': rr.get('m', cr.get('m', 0)),
            'ripopt_status': rr.get('status', 'N/A'),
            'ipopt_status': cr.get('status', 'N/A'),
            'ripopt_obj': rr.get('objective', float('nan')),
            'ipopt_obj': cr.get('objective', float('nan')),
            'obj_diff': od,
            'ripopt_iters': rr.get('iterations', 0),
            'ipopt_iters': cr.get('iterations', 0),
            'ripopt_time': rr.get('solve_time', 0),
            'ipopt_time': cr.get('solve_time', 0),
            'ripopt_solved': r_solved,
            'ipopt_solved': c_solved,
            'both_solved': both,
            'passed': both and not math.isnan(od) and od < 1e-4,
        })

    return comparisons


# ---- Report generation ----

def suite_summary(name, comps):
    """Generate summary stats for a suite."""
    total = len(comps)
    r_solved = sum(1 for c in comps if c['ripopt_solved'])
    i_solved = sum(1 for c in comps if c['ipopt_solved'])
    both = sum(1 for c in comps if c['both_solved'])
    passed = sum(1 for c in comps if c['passed'])

    r_optimal = sum(1 for c in comps if c['ripopt_status'] == 'Optimal')
    r_acceptable = sum(1 for c in comps if c['ripopt_status'] == 'Acceptable')
    i_optimal = sum(1 for c in comps if c['ipopt_status'] == 'Optimal')
    i_acceptable = sum(1 for c in comps if c['ipopt_status'] == 'Acceptable')

    r_only = sum(1 for c in comps if c['ripopt_solved'] and not c['ipopt_solved'])
    i_only = sum(1 for c in comps if c['ipopt_solved'] and not c['ripopt_solved'])

    return {
        'name': name, 'total': total,
        'r_solved': r_solved, 'i_solved': i_solved, 'both': both, 'passed': passed,
        'r_optimal': r_optimal, 'r_acceptable': r_acceptable,
        'i_optimal': i_optimal, 'i_acceptable': i_acceptable,
        'r_only': r_only, 'i_only': i_only,
    }


def speed_stats(comps):
    """Compute speed comparison stats for commonly-solved problems."""
    both_data = [c for c in comps if c['both_solved']
                 and c['ripopt_time'] > 0 and c['ipopt_time'] > 0]
    if not both_data:
        return None

    speedups = [c['ipopt_time'] / c['ripopt_time'] for c in both_data]
    r_times = [c['ripopt_time'] for c in both_data]
    i_times = [c['ipopt_time'] for c in both_data]
    r_iters = [c['ripopt_iters'] for c in both_data]
    i_iters = [c['ipopt_iters'] for c in both_data]

    return {
        'n_problems': len(both_data),
        'geo_mean_speedup': geo_mean(speedups),
        'median_speedup': median(speedups),
        'r_faster_count': sum(1 for s in speedups if s > 1.0),
        'i_faster_count': sum(1 for s in speedups if s < 1.0),
        'r_10x_faster': sum(1 for s in speedups if s > 10.0),
        'r_total_time': sum(r_times),
        'i_total_time': sum(i_times),
        'r_median_time': median(r_times),
        'i_median_time': median(i_times),
        'r_mean_iters': sum(r_iters) / len(r_iters),
        'i_mean_iters': sum(i_iters) / len(i_iters),
        'r_median_iters': median(r_iters),
        'i_median_iters': median(i_iters),
    }


def failure_analysis(comps):
    """Categorize failures by status."""
    r_failures = defaultdict(int)
    i_failures = defaultdict(int)
    for c in comps:
        if not c['ripopt_solved']:
            r_failures[c['ripopt_status']] += 1
        if not c['ipopt_solved']:
            i_failures[c['ipopt_status']] += 1
    return dict(r_failures), dict(i_failures)


def generate_report(suites, output_path, baseline=None):
    """Generate the unified benchmark report."""
    lines = []
    lines.append("# ripopt Benchmark Report")
    lines.append("")
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")

    # Combined summary
    all_comps = []
    for name, comps in suites:
        all_comps.extend(comps)

    combined = suite_summary("Combined", all_comps)

    lines.append("## Executive Summary")
    lines.append("")
    lines.append(f"| Metric | ripopt | Ipopt |")
    lines.append(f"|--------|--------|-------|")
    lines.append(f"| Total solved | **{combined['r_solved']}/{combined['total']}** ({100*combined['r_solved']/max(combined['total'],1):.1f}%) | {combined['i_solved']}/{combined['total']} ({100*combined['i_solved']/max(combined['total'],1):.1f}%) |")
    lines.append(f"| Optimal | {combined['r_optimal']} | {combined['i_optimal']} |")
    lines.append(f"| Acceptable | {combined['r_acceptable']} | {combined['i_acceptable']} |")
    lines.append(f"| Solved exclusively | {combined['r_only']} | {combined['i_only']} |")
    lines.append(f"| Both solved | {combined['both']} | |")
    lines.append(f"| Matching objectives | {combined['passed']}/{max(combined['both'],1)} | |")
    lines.append("")

    # Per-suite summary table
    lines.append("## Per-Suite Summary")
    lines.append("")
    lines.append("| Suite | Problems | ripopt solved | Ipopt solved | ripopt only | Ipopt only | Both solved | Match |")
    lines.append("|-------|----------|--------------|-------------|-------------|------------|------------|-------|")
    for name, comps in suites:
        s = suite_summary(name, comps)
        lines.append(
            f"| {name} | {s['total']} "
            f"| {s['r_solved']} ({100*s['r_solved']/max(s['total'],1):.1f}%) "
            f"| {s['i_solved']} ({100*s['i_solved']/max(s['total'],1):.1f}%) "
            f"| {s['r_only']} | {s['i_only']} | {s['both']} "
            f"| {s['passed']}/{max(s['both'],1)} |"
        )
    lines.append("")

    # Per-suite speed and iteration stats
    for name, comps in suites:
        sp = speed_stats(comps)
        if sp is None:
            continue

        lines.append(f"## {name} Suite — Performance")
        lines.append("")
        lines.append(f"On {sp['n_problems']} commonly-solved problems:")
        lines.append("")
        lines.append("| Metric | ripopt | Ipopt |")
        lines.append("|--------|--------|-------|")
        lines.append(f"| Median time | {fmt_time(sp['r_median_time'])} | {fmt_time(sp['i_median_time'])} |")
        lines.append(f"| Total time | {fmt_time(sp['r_total_time'])} | {fmt_time(sp['i_total_time'])} |")
        lines.append(f"| Mean iterations | {sp['r_mean_iters']:.1f} | {sp['i_mean_iters']:.1f} |")
        lines.append(f"| Median iterations | {sp['r_median_iters']} | {sp['i_median_iters']} |")
        lines.append("")
        lines.append(f"- **Geometric mean speedup**: {sp['geo_mean_speedup']:.1f}x")
        lines.append(f"- **Median speedup**: {sp['median_speedup']:.1f}x")
        lines.append(f"- ripopt faster: {sp['r_faster_count']}/{sp['n_problems']} ({100*sp['r_faster_count']/sp['n_problems']:.0f}%)")
        lines.append(f"- ripopt 10x+ faster: {sp['r_10x_faster']}/{sp['n_problems']}")
        lines.append(f"- Ipopt faster: {sp['i_faster_count']}/{sp['n_problems']}")
        lines.append("")

    # Failure analysis per suite
    lines.append("## Failure Analysis")
    lines.append("")
    for name, comps in suites:
        rf, ifail = failure_analysis(comps)
        if not rf and not ifail:
            continue
        lines.append(f"### {name} Suite")
        lines.append("")
        all_statuses = sorted(set(list(rf.keys()) + list(ifail.keys())))
        lines.append("| Failure Mode | ripopt | Ipopt |")
        lines.append("|-------------|--------|-------|")
        for st in all_statuses:
            lines.append(f"| {st} | {rf.get(st, 0)} | {ifail.get(st, 0)} |")
        lines.append("")

    # Regressions (ripopt fails, ipopt solves)
    regressions = [c for c in all_comps if c['ipopt_solved'] and not c['ripopt_solved']]
    if regressions:
        lines.append("## Regressions (Ipopt solves, ripopt fails)")
        lines.append("")
        lines.append("| Problem | Suite | n | m | ripopt status | Ipopt obj |")
        lines.append("|---------|-------|---|---|--------------|-----------|")
        for c in sorted(regressions, key=lambda c: c['name']):
            io = c['ipopt_obj']
            io_str = f"{io:.6e}" if isinstance(io, (int, float)) and not math.isnan(io) else "N/A"
            lines.append(f"| {c['name']} | {c['suite']} | {c['n']} | {c['m']} | {c['ripopt_status']} | {io_str} |")
        lines.append("")

    # Wins (ripopt solves, ipopt fails)
    wins = [c for c in all_comps if c['ripopt_solved'] and not c['ipopt_solved']]
    if wins:
        lines.append(f"## Wins (ripopt solves, Ipopt fails) — {len(wins)} problems")
        lines.append("")
        lines.append("| Problem | Suite | n | m | Ipopt status | ripopt obj |")
        lines.append("|---------|-------|---|---|-------------|------------|")
        for c in sorted(wins, key=lambda c: c['name']):
            ro = c['ripopt_obj']
            ro_str = f"{ro:.6e}" if isinstance(ro, (int, float)) and not math.isnan(ro) else "N/A"
            lines.append(f"| {c['name']} | {c['suite']} | {c['n']} | {c['m']} | {c['ipopt_status']} | {ro_str} |")
        lines.append("")

    # Acceptable breakdown (problems where ripopt gets Acceptable, not Optimal)
    acceptable = [c for c in all_comps if c['ripopt_status'] == 'Acceptable']
    if acceptable:
        lines.append(f"## Acceptable (not Optimal) — {len(acceptable)} problems")
        lines.append("")
        lines.append("These problems converged within relaxed tolerances but not strict tolerances.")
        lines.append("")
        lines.append("| Problem | Suite | n | m | Ipopt status | ripopt obj | Ipopt obj |")
        lines.append("|---------|-------|---|---|-------------|------------|-----------|")
        for c in sorted(acceptable, key=lambda c: c['name']):
            ro = c['ripopt_obj']
            io = c['ipopt_obj']
            ro_str = f"{ro:.6e}" if isinstance(ro, (int, float)) and not math.isnan(ro) else "N/A"
            io_str = f"{io:.6e}" if isinstance(io, (int, float)) and not math.isnan(io) else "N/A"
            lines.append(f"| {c['name']} | {c['suite']} | {c['n']} | {c['m']} | {c['ipopt_status']} | {ro_str} | {io_str} |")
        lines.append("")

    # Baseline regression detection
    if baseline:
        lines.append("## Regression Detection (vs baseline)")
        lines.append("")
        current_by_name = {c['name']: c for c in all_comps}
        new_failures = []
        new_acceptables = []
        for b in baseline:
            name = b['name']
            if name not in current_by_name:
                continue
            cur = current_by_name[name]
            # Was solved, now fails
            if b['ripopt_solved'] and not cur['ripopt_solved']:
                new_failures.append((name, b['ripopt_status'], cur['ripopt_status']))
            # Was Optimal, now Acceptable
            if b['ripopt_status'] == 'Optimal' and cur['ripopt_status'] == 'Acceptable':
                new_acceptables.append(name)

        if new_failures:
            lines.append(f"### New failures ({len(new_failures)})")
            lines.append("")
            lines.append("| Problem | Was | Now |")
            lines.append("|---------|-----|-----|")
            for name, was, now in new_failures:
                lines.append(f"| {name} | {was} | {now} |")
            lines.append("")

        if new_acceptables:
            lines.append(f"### Degraded to Acceptable ({len(new_acceptables)})")
            lines.append("")
            for name in new_acceptables:
                lines.append(f"- {name}")
            lines.append("")

        if not new_failures and not new_acceptables:
            lines.append("No regressions detected vs baseline.")
            lines.append("")

    # Save machine-readable summary for future regression detection
    summary_data = []
    for c in all_comps:
        summary_data.append({
            'name': c['name'],
            'suite': c['suite'],
            'ripopt_status': c['ripopt_status'],
            'ipopt_status': c['ipopt_status'],
            'ripopt_obj': c['ripopt_obj'] if isinstance(c['ripopt_obj'], (int, float)) and not math.isnan(c['ripopt_obj']) else None,
            'ipopt_obj': c['ipopt_obj'] if isinstance(c['ipopt_obj'], (int, float)) and not math.isnan(c['ipopt_obj']) else None,
            'ripopt_solved': c['ripopt_solved'],
            'ipopt_solved': c['ipopt_solved'],
        })

    lines.append("---")
    lines.append("*Generated by benchmark_report.py*")

    report = '\n'.join(lines)

    with open(output_path, 'w') as f:
        f.write(report)

    # Save baseline JSON for future regression detection
    baseline_path = output_path.replace('.md', '.json')
    with open(baseline_path, 'w') as f:
        json.dump(summary_data, f, indent=2)

    return combined, summary_data


# ---- Main ----

def main():
    output_path = os.path.join(SCRIPT_DIR, 'BENCHMARK_REPORT.md')
    baseline_path = None

    args = sys.argv[1:]
    i = 0
    while i < len(args):
        if args[i] == '--output' and i + 1 < len(args):
            output_path = args[i + 1]
            i += 2
        elif args[i] == '--baseline' and i + 1 < len(args):
            baseline_path = args[i + 1]
            i += 2
        else:
            i += 1

    # Load baseline if provided
    baseline = None
    if baseline_path and os.path.exists(baseline_path):
        with open(baseline_path) as f:
            baseline = json.load(f)
        print(f"Loaded baseline: {baseline_path} ({len(baseline)} problems)")

    # Load all suites
    suites = []

    hs = load_hs_results()
    if hs:
        suites.append(("HS", hs))
        print(f"HS suite: {len(hs)} problems loaded")
    else:
        print("HS suite: no results found (run `make hs-run` first)")

    cutest = load_cutest_results()
    if cutest:
        suites.append(("CUTEst", cutest))
        print(f"CUTEst suite: {len(cutest)} problems loaded")
    else:
        print("CUTEst suite: no results found (run `make cutest-run` first)")

    if not suites:
        print("No benchmark results found. Run `make benchmark` first.")
        sys.exit(1)

    combined, _summary = generate_report(suites, output_path, baseline)

    print(f"\nReport written to {output_path}")
    print(f"Baseline saved to {output_path.replace('.md', '.json')}")
    print(f"\nCombined summary:")
    print(f"  Total: {combined['total']}")
    print(f"  ripopt solved: {combined['r_solved']}/{combined['total']} "
          f"(Optimal: {combined['r_optimal']}, Acceptable: {combined['r_acceptable']})")
    print(f"  Ipopt solved:  {combined['i_solved']}/{combined['total']} "
          f"(Optimal: {combined['i_optimal']}, Acceptable: {combined['i_acceptable']})")
    print(f"  ripopt only:   {combined['r_only']}")
    print(f"  Ipopt only:    {combined['i_only']}")


if __name__ == '__main__':
    main()
