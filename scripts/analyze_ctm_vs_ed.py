#!/usr/bin/env python3
"""Analyze CTM iPEPS trajectory artifacts against the 3x3 PBC ED reference.

Reads the ED CSV at artifacts/m3-systematic/ed-trajectory-3x3-t200-dt002.csv and
all CTM-extended trajectory CSVs under artifacts/m3-systematic-ctm/. For each
(D, cutoff) case, reports max |exact_finite - ED|, max |ctm_density - ED|, and
max |density_simple - ED| up to horizons T in {0.2, 0.5, 1.0, 1.5, 2.0}.

Writes a markdown summary table to stdout and saves it under
artifacts/m3-systematic-ctm/summary.md.
"""
import csv
import glob
import os
import re
import sys


ED_PATH = "artifacts/m3-systematic/ed-trajectory-3x3-t200-dt002.csv"
CTM_DIR = "artifacts/m3-systematic-ctm"
HORIZONS = [0.2, 0.5, 1.0, 1.5, 2.0]


def read_csv_rows(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def parse_case_from_filename(path):
    m = re.search(r"D(\d+)-cut([0-9eE.+-]+)\.csv$", path)
    if not m:
        return None
    D = int(m.group(1))
    cutoff = float(m.group(2))
    return D, cutoff


def float_or_none(s):
    if s is None or s == "":
        return None
    try:
        v = float(s)
    except ValueError:
        return None
    return v


def max_abs_up_to(values, T):
    return max((abs(v) for v in values if v[0] is None or v[0] <= T + 1e-9), default=None)


def main():
    if not os.path.isfile(ED_PATH):
        sys.exit(f"missing ED reference: {ED_PATH}")
    ed_rows = read_csv_rows(ED_PATH)
    ed_density = {float(r["time"]): float(r["excitation_density"]) for r in ed_rows}

    files = sorted(glob.glob(os.path.join(CTM_DIR, "*.csv")))
    if not files:
        sys.exit(f"no CTM trajectory CSVs in {CTM_DIR}")

    rows = []
    for f in files:
        case = parse_case_from_filename(f)
        if case is None:
            continue
        D, cutoff = case
        traj = read_csv_rows(f)
        errs_ctm = {T: [] for T in HORIZONS}
        errs_exact = {T: [] for T in HORIZONS}
        errs_simple = {T: [] for T in HORIZONS}
        ctm_seconds = []
        max_ctm_residual = 0.0
        all_accepted = True
        for r in traj:
            t = float(r["time"])
            if t not in ed_density:
                continue
            ed_v = ed_density[t]
            ctm_v = float_or_none(r.get("ctm_density"))
            exa_v = float_or_none(r.get("exact_finite_density"))
            sim_v = float_or_none(r.get("density_simple"))
            res = float_or_none(r.get("ctm_residual"))
            acc = r.get("ctm_accepted", "")
            secs = float_or_none(r.get("ctm_seconds"))
            if secs is not None:
                ctm_seconds.append(secs)
            if res is not None and res > max_ctm_residual:
                max_ctm_residual = res
            if acc not in ("", "true", "True", None):
                all_accepted = False
            for T in HORIZONS:
                if t <= T + 1e-9:
                    if ctm_v is not None:
                        errs_ctm[T].append(abs(ctm_v - ed_v))
                    if exa_v is not None:
                        errs_exact[T].append(abs(exa_v - ed_v))
                    if sim_v is not None:
                        errs_simple[T].append(abs(sim_v - ed_v))
        row = {
            "D": D,
            "cutoff": cutoff,
            "ctm_total_seconds": sum(ctm_seconds),
            "max_ctm_residual": max_ctm_residual,
            "all_accepted": all_accepted,
        }
        for T in HORIZONS:
            row[f"err_ctm_T{T}"] = max(errs_ctm[T], default=None)
            row[f"err_exact_T{T}"] = max(errs_exact[T], default=None)
            row[f"err_simple_T{T}"] = max(errs_simple[T], default=None)
        rows.append(row)

    rows.sort(key=lambda r: (r["D"], r["cutoff"]))

    def fmt(x):
        return "    -    " if x is None else f"{x:.2e}"

    out = []
    out.append("# CTM vs ED on 3x3 PBC PXP (all-down, dt=0.02, T=2.0)")
    out.append("")
    out.append("Three observables per case at increasing time horizons T:")
    out.append("- err_ctm: max |ctm_density - ed_density|")
    out.append("- err_exact: max |exact_finite_density - ed_density|")
    out.append("- err_simple: max |density_simple - ed_density|")
    out.append("")
    for label in ("err_ctm", "err_exact", "err_simple"):
        out.append(f"## {label}")
        out.append("")
        header = "| D | cutoff | " + " | ".join(f"T={T}" for T in HORIZONS) + " |"
        sep = "|---|---|" + "|".join(["---:"] * len(HORIZONS)) + "|"
        out.append(header)
        out.append(sep)
        for r in rows:
            cells = " | ".join(fmt(r[f"{label}_T{T}"]) for T in HORIZONS)
            out.append(f"| {r['D']} | {r['cutoff']:.0e} | {cells} |")
        out.append("")

    out.append("## CTM diagnostics")
    out.append("")
    out.append("| D | cutoff | total CTM seconds | max ctm_residual | all accepted |")
    out.append("|---|---|---:|---:|---|")
    for r in rows:
        out.append(
            f"| {r['D']} | {r['cutoff']:.0e} | {r['ctm_total_seconds']:.1f} | "
            f"{r['max_ctm_residual']:.2e} | {r['all_accepted']} |"
        )
    out.append("")

    md = "\n".join(out)
    print(md)
    summary_path = os.path.join(CTM_DIR, "summary.md")
    with open(summary_path, "w") as f:
        f.write(md + "\n")
    print(f"\nWrote {summary_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
