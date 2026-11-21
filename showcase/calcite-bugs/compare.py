"""Run the original calcite model and every fix, and plot them side by side.

    python compare.py [--phreeqc PATH] [--database PATH]

original/calcite.phr is models/calcite/calcite.phr with the placeholders
filled in with the smoke-test values (the median of every parameter in
models/calcite/config.cfg). The variants are built from it with the edits in
fixes.py: each fix on its own, and all fixes together. See README.md.

Writes, next to this script:
    variants/<name>/calcite.phr     the model of each variant ("# FIX n:" marks the edits)
    output/<name>.tsv               PHREEQC's selected output of each run
    figures/<fix>.png               one fix vs the original
    figures/all_fixes.png           the original, every fix and all fixes together
"""

import argparse
import configparser
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from fixes import FIXES, apply

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ORIGINAL = HERE / "original" / "calcite.phr"

# name -> (column title, fixes applied)
VARIANTS = {"original": ("Original", [])}
VARIANTS.update({name: (fix["title"], [name]) for name, fix in FIXES.items()})
VARIANTS["all-fixes"] = ("All fixes", list(FIXES))

# One row per panel: (title, columns, scale, ylabel). "reacted:" columns are
# kinetic reactants, plotted as the moles reacted so far (-m minus KIN()).
PANELS = [
    ("pH", ["pH"], 1, "pH"),
    ("Gas and solution", ["H2_g_mol", "CO2_g_mol", "SO4_mol", "Formate", "Acetate", "Ca_mol"], 1000, "Amount (mmol)"),
    ("Biomass", ["X_SRB_mol", "X_IRB_mol"], 1, "Moles"),
    ("Ca, Fe, HS-", ["Ca_mol", "Fe_mol", "HS-_mol"], 1000, "Amount (mmol)"),
    ("Sulfate reduction (reacted)", ["reacted:SR_H2", "reacted:SR_AC", "reacted:SR_FOR"], 1e6, "Reacted (µmol)"),
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


def run(path, phreeqc, database):
    with tempfile.TemporaryDirectory() as tmp:
        pqo = Path(tmp, "model.pqo")
        res = subprocess.run(
            [phreeqc, str(path), str(pqo), database],
            capture_output=True, text=True, cwd=tmp, check=False,
        )
        out = pqo.read_text(encoding="utf-8") if pqo.exists() else ""
        errors = [ln for ln in out.splitlines() if ln.startswith("ERROR") and ln.strip() != "ERROR:"]
        if res.returncode != 0 or errors:
            sys.exit(f"{path}: PHREEQC failed: {' | '.join(errors[:3]) or res.stderr.strip()}")
        sel = Path(tmp, "selected_output_1.sel").read_text(encoding="utf-8")
        version = re.search(r"\*\s*(PHREEQC[\w.\-]*)\s*\*", res.stdout + res.stderr + out)
    return sel, version.group(1) if version else "PHREEQC"


def parse(sel):
    lines = [ln for ln in sel.splitlines() if ln.strip()]
    col = {h: i for i, h in enumerate(lines[0].split())}
    data = np.array([[float(v) for v in ln.split()] for ln in lines[1:]])
    data = data[np.argmax(data[:, col["time_d"]] > 0):]  # drop the t = 0 initial-solution row
    return col, data


def series(result, c):
    col, data, m0 = result
    if c.startswith("reacted:"):
        r = c.split(":", 1)[1]
        return f"{r} reacted", m0[r] - data[:, col[r]]
    return c, data[:, col[c]]


def figure(results, names, title, path):
    """One column per variant, one row per panel, the same y axis along a row."""
    ncol = len(names)
    fig, axes = plt.subplots(len(PANELS), ncol, figsize=(5.2 * ncol + 1, 3.4 * len(PANELS)),
                             sharex=True, squeeze=False)
    fig.suptitle(title, fontsize=15, fontweight="bold")
    for row, (panel, columns, scale, ylabel) in enumerate(PANELS):
        for j, name in enumerate(names):
            ax = axes[row, j]
            t = results[name][1][:, results[name][0]["time_d"]]
            for c in columns:
                label, y = series(results[name], c)
                ax.plot(t, y * scale, label=label.replace(" reacted", ""))
            ax.set_title(f"{VARIANTS[name][0]}\n{panel}" if row == 0 else panel, fontsize=10)
            ax.grid(True, linestyle="--", alpha=0.6)
            if j == 0:
                ax.set_ylabel(ylabel)
            if len(columns) > 1 and j == ncol - 1:
                ax.legend(loc="best", fontsize="x-small", ncol=2)
        lo = min(a.get_ylim()[0] for a in axes[row])
        hi = max(a.get_ylim()[1] for a in axes[row])
        for a in axes[row]:
            a.set_ylim(lo, hi)
    for a in axes[-1]:
        a.set_xlabel("Time (days)")
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(path, dpi=110)
    plt.close(fig)
    print(f"Figure: {path.relative_to(HERE)}")


def table(results, names):
    width = 15
    print(f"\n{'end of run':<20}" + "".join(f"{n:>{width}}" for n in names))
    seen = set()
    for _, columns, _, _ in PANELS:
        for c in columns:
            label = series(results[names[0]], c)[0]
            if label in seen:
                continue
            seen.add(label)
            print(f"{label:<20}" + "".join(f"{series(results[n], c)[1][-1]:>{width}.4e}" for n in names))


def main():
    exe, db = default_phreeqc()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phreeqc", default=exe, help=f"PHREEQC executable (default from config.cfg: {exe})")
    ap.add_argument("--database", default=db, help=f"phreeqc.dat (default from config.cfg: {db})")
    args = ap.parse_args()

    original = ORIGINAL.read_text(encoding="utf-8")
    paths = {"original": ORIGINAL}
    for name, (_, fixes) in VARIANTS.items():
        if fixes:
            paths[name] = HERE / "variants" / name / "calcite.phr"
            paths[name].parent.mkdir(parents=True, exist_ok=True)
            paths[name].write_text(apply(original, fixes), encoding="utf-8")

    print(f"Running {len(paths)} models ...", flush=True)
    with ThreadPoolExecutor() as pool:
        runs = dict(zip(paths, pool.map(lambda p: run(p, args.phreeqc, args.database), paths.values())))

    (HERE / "output").mkdir(exist_ok=True)
    results = {}
    for name, (sel, version) in runs.items():
        (HERE / "output" / f"{name}.tsv").write_text(sel, encoding="utf-8")
        col, data = parse(sel)
        results[name] = (col, data, initial_moles(paths[name].read_text(encoding="utf-8")))

    figs = HERE / "figures"
    figs.mkdir(exist_ok=True)
    for name, fix in FIXES.items():
        figure(results, ["original", name],
               f"Calcite model: original vs {fix['title']}\n{fix['summary']} ({version})",
               figs / f"{name}.png")
    figure(results, list(VARIANTS), f"Calcite model: original, each fix alone and all fixes ({version})",
           figs / "all_fixes.png")
    table(results, list(VARIANTS))


if __name__ == "__main__":
    main()
