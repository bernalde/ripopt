#!/usr/bin/env python3
"""Generate a markdown comparison report for the Mittelmann ampl-nlp benchmark.

Reads results/ripopt_<version>.json and results/ipopt_<version>.json (if present)
and writes a Mittelmann-style table:

  problem | n | m | ripopt status | ripopt time | ipopt status | ipopt time
  ...
  totals  |   |   | <count>       | SGM         | <count>      | SGM

SGM = scaled shifted geometric mean, shift=10s, matching Mittelmann's convention.
Failures (TIMEOUT, ERROR) are scored at the time limit.
"""
import json
import math
import sys
from pathlib import Path


def sgm(times, shift=10.0):
    """Scaled shifted geometric mean."""
    if not times:
        return float("nan")
    return math.exp(sum(math.log(t + shift) for t in times) / len(times)) - shift


def load(version, solver):
    p = Path(f"results/{solver}_{version}.json")
    if not p.exists():
        return {}
    with open(p) as f:
        rows = json.load(f)
    return {r["problem"]: r for r in rows}


def main():
    if len(sys.argv) != 3:
        print("usage: make_report.py <version> <out.md>", file=sys.stderr)
        sys.exit(1)
    version, out_path = sys.argv[1], sys.argv[2]

    rip = load(version, "ripopt")
    ipo = load(version, "ipopt")
    problems = sorted(set(rip) | set(ipo))

    lines = [
        f"# Mittelmann ampl-nlp benchmark — ripopt {version}",
        "",
        "Source: https://plato.asu.edu/ftp/ampl-nlp.html",
        "",
        "Status codes: `Optimal` = solver returned success; `TIMEOUT` = exceeded "
        "wall-clock limit; `ERROR` = nonzero exit. Times are wall-clock seconds.",
        "",
        "| problem | ripopt status | ripopt time | ripopt iter | ipopt status | ipopt time | ipopt iter |",
        "|---------|---------------|-------------|-------------|--------------|------------|------------|",
    ]

    rip_times, ipo_times = [], []
    rip_solved = ipo_solved = 0
    for p in problems:
        r = rip.get(p)
        i = ipo.get(p)
        rs = r["status"] if r else "—"
        is_ = i["status"] if i else "—"
        rt = r["elapsed"] if r else None
        it = i["elapsed"] if i else None
        ri = r["iterations"] if r else None
        ii = i["iterations"] if i else None
        if r and rs == "OK":
            rip_solved += 1
            rip_times.append(rt or 0.0)
        if i and is_ == "OK":
            ipo_solved += 1
            ipo_times.append(it or 0.0)
        lines.append(
            f"| {p} | {rs} | {rt or '—'} | {ri or '—'} | {is_} | {it or '—'} | {ii or '—'} |"
        )

    lines.append("")
    lines.append(f"**Solved:** ripopt {rip_solved}/{len(problems)}  |  ipopt {ipo_solved}/{len(problems)}")
    lines.append(
        f"**SGM (shift=10s, solved only):** ripopt {sgm(rip_times):.2f}  |  ipopt {sgm(ipo_times):.2f}"
    )

    Path(out_path).write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
