#!/usr/bin/env python3
"""
Compare ripopt and cyipopt results on HS test problems.
Generates HS_VALIDATION_REPORT.md.
"""

import json
import os
import sys
import math
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def load_results(path):
    """Load JSON results and index by problem number."""
    with open(path) as f:
        data = json.load(f)
    return {r['number']: r for r in data}


def is_solved(r):
    return r['status'] in ('Optimal', 'Acceptable')


def obj_diff(ripopt_r, cyipopt_r):
    """Relative objective difference between solvers.

    Uses max(|f_r|, |f_c|, 1.0) as denominator to handle near-zero objectives.
    """
    ro = ripopt_r.get('objective')
    co = cyipopt_r.get('objective')
    if ro is None or co is None:
        return float('nan')
    if math.isnan(ro) or math.isnan(co):
        return float('nan')
    denom = max(abs(co), abs(ro), 1.0)
    return abs(ro - co) / denom


def x_diff(ripopt_r, cyipopt_r):
    """Max absolute difference in solution vectors."""
    rx = ripopt_r.get('x', [])
    cx = cyipopt_r.get('x', [])
    if not rx or not cx or len(rx) != len(cx):
        return float('nan')
    try:
        return max(abs(a - b) for a, b in zip(rx, cx) if a is not None and b is not None)
    except (TypeError, ValueError):
        return float('nan')


def fmt_time(t):
    """Format a time value in seconds to a human-readable string."""
    if t >= 1.0:
        return f"{t:.2f}s"
    elif t >= 0.001:
        return f"{t*1000:.1f}ms"
    else:
        return f"{t*1e6:.0f}us"


def classify_problem(r):
    """Classify by constraint type."""
    n = r.get('n', 0)
    m = r.get('m', 0)
    if m == 0:
        return "unconstrained"
    else:
        return "constrained"


def main():
    ripopt_path = os.path.join(SCRIPT_DIR, 'hs_ripopt_results.json')
    cyipopt_path = os.path.join(SCRIPT_DIR, 'hs_cyipopt_results.json')
    report_path = os.path.join(SCRIPT_DIR, 'HS_VALIDATION_REPORT.md')

    ripopt = load_results(ripopt_path)
    cyipopt = load_results(cyipopt_path)

    all_numbers = sorted(set(ripopt.keys()) | set(cyipopt.keys()))

    # Compute comparison metrics
    comparisons = []
    for num in all_numbers:
        rr = ripopt.get(num, {})
        cr = cyipopt.get(num, {})

        r_solved = is_solved(rr) if rr else False
        c_solved = is_solved(cr) if cr else False
        both_solved = r_solved and c_solved

        od = obj_diff(rr, cr) if both_solved else float('nan')
        xd = x_diff(rr, cr) if both_solved else float('nan')

        # Pass criteria: both solve and obj difference < 1e-4
        if both_solved and not math.isnan(od):
            passed = od < 1e-4
        else:
            passed = False

        comparisons.append({
            'number': num,
            'n': rr.get('n', cr.get('n', 0)),
            'm': rr.get('m', cr.get('m', 0)),
            'ripopt_status': rr.get('status', 'N/A'),
            'cyipopt_status': cr.get('status', 'N/A'),
            'ripopt_obj': rr.get('objective', float('nan')),
            'cyipopt_obj': cr.get('objective', float('nan')),
            'obj_diff': od,
            'x_diff': xd,
            'ripopt_iters': rr.get('iterations', 0),
            'cyipopt_iters': cr.get('iterations', 0),
            'ripopt_time': rr.get('solve_time', 0),
            'cyipopt_time': cr.get('solve_time', 0),
            'known_fopt': rr.get('known_fopt', cr.get('known_fopt', float('nan'))),
            'passed': passed,
            'ripopt_solved': r_solved,
            'cyipopt_solved': c_solved,
            'both_solved': both_solved,
            'category': classify_problem(rr if rr else cr),
        })

    # Statistics
    total = len(comparisons)
    ripopt_solved = sum(1 for c in comparisons if c['ripopt_solved'])
    cyipopt_solved = sum(1 for c in comparisons if c['cyipopt_solved'])
    both_solved = sum(1 for c in comparisons if c['both_solved'])
    passed = sum(1 for c in comparisons if c['passed'])

    matched_diffs = [c['obj_diff'] for c in comparisons if c['both_solved'] and not math.isnan(c['obj_diff'])]
    matched_xdiffs = [c['x_diff'] for c in comparisons if c['both_solved'] and not math.isnan(c['x_diff'])]

    if matched_diffs:
        mean_obj_diff = sum(matched_diffs) / len(matched_diffs)
        max_obj_diff = max(matched_diffs)
        median_obj_diff = sorted(matched_diffs)[len(matched_diffs) // 2]
    else:
        mean_obj_diff = max_obj_diff = median_obj_diff = float('nan')

    if matched_xdiffs:
        mean_x_diff = sum(matched_xdiffs) / len(matched_xdiffs)
        max_x_diff = max(matched_xdiffs)
        median_x_diff = sorted(matched_xdiffs)[len(matched_xdiffs) // 2]
    else:
        mean_x_diff = max_x_diff = median_x_diff = float('nan')

    # Category breakdown
    categories = defaultdict(lambda: {'total': 0, 'ripopt': 0, 'cyipopt': 0, 'both': 0, 'passed': 0})
    for c in comparisons:
        cat = c['category']
        categories[cat]['total'] += 1
        if c['ripopt_solved']:
            categories[cat]['ripopt'] += 1
        if c['cyipopt_solved']:
            categories[cat]['cyipopt'] += 1
        if c['both_solved']:
            categories[cat]['both'] += 1
        if c['passed']:
            categories[cat]['passed'] += 1

    # Only ripopt fails
    ripopt_only_fails = [c for c in comparisons if not c['ripopt_solved'] and c['cyipopt_solved']]
    # Only cyipopt fails
    cyipopt_only_fails = [c for c in comparisons if c['ripopt_solved'] and not c['cyipopt_solved']]
    # Both fail
    both_fail = [c for c in comparisons if not c['ripopt_solved'] and not c['cyipopt_solved']]

    # Generate report
    lines = []
    lines.append("# HS Test Suite Validation Report")
    lines.append("")
    lines.append("Comparison of ripopt vs cyipopt (C++ Ipopt) on the Hock-Schittkowski test suite.")
    lines.append("")
    lines.append("## Executive Summary")
    lines.append("")
    lines.append(f"- **Total problems**: {total}")
    lines.append(f"- **ripopt solved**: {ripopt_solved}/{total} ({100*ripopt_solved/total:.1f}%)")
    lines.append(f"- **cyipopt solved**: {cyipopt_solved}/{total} ({100*cyipopt_solved/total:.1f}%)")
    lines.append(f"- **Both solved**: {both_solved}/{total} ({100*both_solved/total:.1f}%)")
    lines.append(f"- **Matching solutions** (rel obj diff < 1e-4): {passed}/{both_solved} ({100*passed/max(both_solved,1):.1f}%)")
    lines.append("")

    lines.append("## Accuracy Statistics (where both solve)")
    lines.append("")
    lines.append(f"| Metric | Objective Rel Diff | Solution Max Diff |")
    lines.append(f"|--------|-------------------|-------------------|")
    lines.append(f"| Mean   | {mean_obj_diff:.2e} | {mean_x_diff:.2e} |")
    lines.append(f"| Median | {median_obj_diff:.2e} | {median_x_diff:.2e} |")
    lines.append(f"| Max    | {max_obj_diff:.2e} | {max_x_diff:.2e} |")
    lines.append("")

    lines.append("## Category Breakdown")
    lines.append("")
    lines.append("| Category | Total | ripopt | cyipopt | Both | Match |")
    lines.append("|----------|-------|--------|---------|------|-------|")
    for cat in sorted(categories.keys()):
        d = categories[cat]
        lines.append(f"| {cat} | {d['total']} | {d['ripopt']} | {d['cyipopt']} | {d['both']} | {d['passed']} |")
    lines.append("")

    lines.append("## Detailed Results")
    lines.append("")
    lines.append("| TP# | n | m | ripopt | cyipopt | Obj Diff | r_iter | c_iter | r_time | c_time | Speedup | Status |")
    lines.append("|-----|---|---|--------|---------|----------|--------|--------|--------|--------|---------|--------|")
    for c in comparisons:
        num = c['number']
        n = c['n']
        m = c['m']
        rs = c['ripopt_status'][:12]
        cs = c['cyipopt_status'][:12]
        od = f"{c['obj_diff']:.2e}" if not math.isnan(c['obj_diff']) else "N/A"
        ri = c['ripopt_iters']
        ci = c['cyipopt_iters']
        rt = c['ripopt_time']
        ct = c['cyipopt_time']
        rt_str = fmt_time(rt) if rt > 0 else "N/A"
        ct_str = fmt_time(ct) if ct > 0 else "N/A"
        if rt > 0 and ct > 0:
            speedup = ct / rt
            sp_str = f"{speedup:.1f}x"
        else:
            sp_str = "N/A"
        status = "PASS" if c['passed'] else ("BOTH_FAIL" if not c['ripopt_solved'] and not c['cyipopt_solved'] else ("ripopt_FAIL" if not c['ripopt_solved'] else "MISMATCH"))
        lines.append(f"| {num:03d} | {n} | {m} | {rs} | {cs} | {od} | {ri} | {ci} | {rt_str} | {ct_str} | {sp_str} | {status} |")
    lines.append("")

    lines.append("## Performance Comparison (where both solve)")
    lines.append("")

    if both_solved > 0:
        both_data = [c for c in comparisons if c['both_solved']]
        r_iters = [c['ripopt_iters'] for c in both_data]
        c_iters = [c['cyipopt_iters'] for c in both_data]
        r_times = [c['ripopt_time'] for c in both_data if c['ripopt_time'] > 0]
        c_times = [c['cyipopt_time'] for c in both_data if c['cyipopt_time'] > 0]

        lines.append("### Iteration Comparison")
        lines.append("")
        lines.append(f"| Metric | ripopt | cyipopt |")
        lines.append(f"|--------|--------|---------|")
        lines.append(f"| Mean   | {sum(r_iters)/len(r_iters):.1f} | {sum(c_iters)/len(c_iters):.1f} |")
        lines.append(f"| Median | {sorted(r_iters)[len(r_iters)//2]} | {sorted(c_iters)[len(c_iters)//2]} |")
        lines.append(f"| Max    | {max(r_iters)} | {max(c_iters)} |")
        lines.append(f"| Total  | {sum(r_iters)} | {sum(c_iters)} |")
        lines.append("")

        # Count how often each solver uses fewer iterations
        r_fewer = sum(1 for r, c in zip(r_iters, c_iters) if r < c)
        c_fewer = sum(1 for r, c in zip(r_iters, c_iters) if c < r)
        tied = sum(1 for r, c in zip(r_iters, c_iters) if r == c)
        lines.append(f"- ripopt uses fewer iterations: {r_fewer}/{len(r_iters)} problems")
        lines.append(f"- cyipopt uses fewer iterations: {c_fewer}/{len(c_iters)} problems")
        lines.append(f"- Same iteration count: {tied}/{len(r_iters)} problems")
        lines.append("")

        if r_times and c_times:
            lines.append("### Timing Comparison")
            lines.append("")
            r_total = sum(r_times)
            c_total = sum(c_times)
            lines.append(f"| Metric | ripopt | cyipopt |")
            lines.append(f"|--------|--------|---------|")
            lines.append(f"| Mean   | {fmt_time(sum(r_times)/len(r_times))} | {fmt_time(sum(c_times)/len(c_times))} |")
            lines.append(f"| Median | {fmt_time(sorted(r_times)[len(r_times)//2])} | {fmt_time(sorted(c_times)[len(c_times)//2])} |")
            lines.append(f"| Max    | {fmt_time(max(r_times))} | {fmt_time(max(c_times))} |")
            lines.append(f"| Total  | {fmt_time(r_total)} | {fmt_time(c_total)} |")
            lines.append("")

            # Geometric mean speedup (cyipopt_time / ripopt_time)
            speedups = [c['cyipopt_time'] / c['ripopt_time']
                        for c in both_data
                        if c['ripopt_time'] > 0 and c['cyipopt_time'] > 0]
            if speedups:
                geo_mean = math.exp(sum(math.log(s) for s in speedups) / len(speedups))
                r_faster = sum(1 for s in speedups if s > 1.0)
                c_faster = sum(1 for s in speedups if s < 1.0)
                lines.append(f"- Geometric mean speedup (cyipopt_time/ripopt_time): **{geo_mean:.2f}x**")
                lines.append(f"  - \\>1 means ripopt is faster, <1 means cyipopt is faster")
                lines.append(f"- ripopt faster: {r_faster}/{len(speedups)} problems")
                lines.append(f"- cyipopt faster: {c_faster}/{len(speedups)} problems")
                lines.append(f"- Overall speedup (total time): {c_total/r_total:.2f}x")
                lines.append("")

    lines.append("")

    lines.append("## Failure Analysis")
    lines.append("")

    if ripopt_only_fails:
        lines.append(f"### Problems where only ripopt fails ({len(ripopt_only_fails)})")
        lines.append("")
        lines.append("| TP# | n | m | ripopt status | cyipopt obj |")
        lines.append("|-----|---|---|---------------|-------------|")
        for c in ripopt_only_fails:
            lines.append(f"| {c['number']:03d} | {c['n']} | {c['m']} | {c['ripopt_status']} | {c['cyipopt_obj']:.6e} |")
        lines.append("")

    if cyipopt_only_fails:
        lines.append(f"### Problems where only cyipopt fails ({len(cyipopt_only_fails)})")
        lines.append("")
        lines.append("| TP# | n | m | cyipopt status | ripopt obj |")
        lines.append("|-----|---|---|----------------|------------|")
        for c in cyipopt_only_fails:
            ro = c['ripopt_obj']
            ro_str = f"{ro:.6e}" if ro is not None and not math.isnan(ro) else "N/A"
            lines.append(f"| {c['number']:03d} | {c['n']} | {c['m']} | {c['cyipopt_status']} | {ro_str} |")
        lines.append("")

    if both_fail:
        lines.append(f"### Problems where both fail ({len(both_fail)})")
        lines.append("")
        lines.append("| TP# | n | m | ripopt status | cyipopt status |")
        lines.append("|-----|---|---|---------------|----------------|")
        for c in both_fail:
            lines.append(f"| {c['number']:03d} | {c['n']} | {c['m']} | {c['ripopt_status']} | {c['cyipopt_status']} |")
        lines.append("")

    # Mismatches where both solve but disagree
    mismatches = [c for c in comparisons if c['both_solved'] and not c['passed']]
    if mismatches:
        lines.append(f"### Objective mismatches (both solve but differ > 1e-4) ({len(mismatches)})")
        lines.append("")
        lines.append("| TP# | ripopt obj | cyipopt obj | Rel Diff |")
        lines.append("|-----|-----------|-------------|----------|")
        def _fmt(v):
            if v is None or (isinstance(v, float) and math.isnan(v)):
                return "N/A"
            return f"{v:.6e}"
        for c in mismatches:
            lines.append(f"| {c['number']:03d} | {_fmt(c['ripopt_obj'])} | {_fmt(c['cyipopt_obj'])} | {c['obj_diff']:.2e} |")
        lines.append("")

    lines.append("---")
    lines.append(f"*Generated by hs_suite/compare.py*")

    report = '\n'.join(lines)

    with open(report_path, 'w') as f:
        f.write(report)

    print(f"Report written to {report_path}")
    print(f"\nSummary:")
    print(f"  Total: {total}")
    print(f"  ripopt solved: {ripopt_solved}/{total}")
    print(f"  cyipopt solved: {cyipopt_solved}/{total}")
    print(f"  Both solved: {both_solved}/{total}")
    print(f"  Matching (rel diff < 1e-4): {passed}/{both_solved}")


if __name__ == '__main__':
    main()
