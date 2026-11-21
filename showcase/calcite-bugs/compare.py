"""Run the original and the fixed calcite model and plot them side by side.

    python compare.py [--phreeqc PATH] [--database PATH]

Both models are models/calcite/calcite.phr with the placeholders filled in
with the smoke-test values (the median of every parameter in
models/calcite/config.cfg). fixed/calcite.phr also has the bug fixes, each
marked with "# FIX:". See README.md.

Writes, next to this script:
    output/original.tsv, output/fixed.tsv   PHREEQC's selected output
    calcite_comparison.png                  the figure
"""

import argparse
import configparser
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
MODELS = {"original": HERE / "original" / "calcite.phr", "fixed": HERE / "fixed" / "calcite.phr"}

# One row per panel: (title, columns, scale, ylabel). "reacted:" columns are
# kinetic reactants, plotted as the moles reacted so far (-m minus KIN()).
PANELS = [
    ("pH", ["pH"], 1, "pH"),
    ("Gas and solution", ["H2_g_mol", "CO2_g_mol", "SO4_mol", "Formate", "Acetate", "Ca_mol"], 1000, "Amount (mmol)"),
    ("Biomass", ["X_SRB_mol", "X_IRB_mol"], 1, "Moles"),
    ("Ca, Fe, HS-", ["Ca_mol", "Fe_mol", "HS-_mol"], 1000, "Amount (mmol)"),
    ("Sulfate reduction (reacted)", ["reacted:SR_H2", "reacted:SR_AC", "reacted:SR_FOR"], 1000, "Reacted (mmol)"),
    ("Iron reduction (reacted)", ["reacted:IR_H2", "reacted:IR_AC", "reacted:IR_FOR"], 1e6, "Reacted (µmol)"),
    ("Gas volume", ["Gas_VOL"], 1, "Gas volume (L)"),
]


def default_phreeqc():
    cp = configparser.ConfigParser(inline_comment_prefixes=("#",))
    cp.read(REPO / "models" / "calcite" / "config.cfg", encoding="utf-8")
    return cp["phreeqc"]["executable"], cp["phreeqc"]["database"]


def initial_moles(template):
    """-m of every kinetic reactant in the KINETICS block."""
    block = template.split("KINETICS 1", 1)[1].split("RATES", 1)[0]
    m0, name = {}, None
    for line in block.splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        if not line.startswith("-"):
            name = line.split()[0]
        elif line.startswith("-m ") and name:
            m0[name] = float(line.split()[1])
    return m0


def run(name, phreeqc, database):
    with tempfile.TemporaryDirectory() as tmp:
        pqo = Path(tmp, "model.pqo")
        res = subprocess.run(
            [phreeqc, str(MODELS[name]), str(pqo), database],
            capture_output=True, text=True, cwd=tmp, check=False,
        )
        out = pqo.read_text(encoding="utf-8") if pqo.exists() else ""
        errors = [ln for ln in out.splitlines() if ln.startswith("ERROR") and ln.strip() != "ERROR:"]
        if res.returncode != 0 or errors:
            sys.exit(f"{name}: PHREEQC failed: {' | '.join(errors[:3]) or res.stderr.strip()}")
        sel = Path(tmp, "selected_output_1.sel").read_text(encoding="utf-8")
        version = re.search(r"\*\s*(PHREEQC[\w.\-]*)\s*\*", res.stdout + res.stderr + out)
    return sel, version.group(1) if version else "PHREEQC"


def parse(sel):
    lines = [ln for ln in sel.splitlines() if ln.strip()]
    headers = lines[0].split()
    data = np.array([[float(v) for v in ln.split()] for ln in lines[1:]])
    col = {h: i for i, h in enumerate(headers)}
    data = data[np.argmax(data[:, col["time_d"]] > 0):]  # drop the t = 0 initial-solution row
    return col, data


def series(col, data, m0, c):
    if c.startswith("reacted:"):
        r = c.split(":", 1)[1]
        return r, m0[r] - data[:, col[r]]
    return c, data[:, col[c]]


def main():
    exe, db = default_phreeqc()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phreeqc", default=exe, help=f"PHREEQC executable (default from config.cfg: {exe})")
    ap.add_argument("--database", default=db, help=f"phreeqc.dat (default from config.cfg: {db})")
    args = ap.parse_args()

    (HERE / "output").mkdir(exist_ok=True)
    results = {}
    for name in MODELS:
        print(f"Running {name} ...", flush=True)
        sel, version = run(name, args.phreeqc, args.database)
        (HERE / "output" / f"{name}.tsv").write_text(sel, encoding="utf-8")
        col, data = parse(sel)
        m0 = initial_moles(MODELS[name].read_text(encoding="utf-8"))
        results[name] = (col, data, m0)

    fig, axes = plt.subplots(len(PANELS), 2, figsize=(14, 3.6 * len(PANELS)), sharex=True, squeeze=False)
    fig.suptitle(f"Calcite model, default values: original vs fixed ({version})", fontsize=14, fontweight="bold")
    for row, (title, columns, scale, ylabel) in enumerate(PANELS):
        for j, name in enumerate(MODELS):
            col, data, m0 = results[name]
            ax = axes[row, j]
            t = data[:, col["time_d"]]
            for c in columns:
                label, y = series(col, data, m0, c)
                ax.plot(t, y * scale, label=label)
            ax.set_title(f"{name.capitalize()}: {title}")
            ax.set_ylabel(ylabel)
            ax.grid(True, linestyle="--", alpha=0.6)
            if len(columns) > 1:
                ax.legend(loc="best", fontsize="small", ncol=2)
        # Same y axis in both columns of a row, so the difference is visible.
        lo = min(a.get_ylim()[0] for a in axes[row]); hi = max(a.get_ylim()[1] for a in axes[row])
        for a in axes[row]:
            a.set_ylim(lo, hi)
    for a in axes[-1]:
        a.set_xlabel("Time (days)")
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    png = HERE / "calcite_comparison.png"
    fig.savefig(png, dpi=130)
    print(f"Figure: {png}")

    # Final values of every plotted quantity.
    print(f"\n{'quantity (end of run)':<24}{'original':>14}{'fixed':>14}")
    seen = set()
    for _, columns, _, _ in PANELS:
        for c in columns:
            vals = [series(*results[n][:2], results[n][2], c) for n in MODELS]
            label = vals[0][0] + (" reacted" if c.startswith("reacted:") else "")
            if label in seen:
                continue
            seen.add(label)
            print(f"{label:<24}{vals[0][1][-1]:>14.4e}{vals[1][1][-1]:>14.4e}")


if __name__ == "__main__":
    main()
