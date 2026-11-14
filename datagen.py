"""
Synthetic data generation for PHREEQC models.

Runs every enabled model in the models/ folder next to this script:

    python datagen.py

Each model has its own folder, named after the model:

    models/<name>/
        <name>.phr                PHREEQC template
        config.cfg                its settings (enabled = false skips the model)
        runs/<timestamp>/         one folder per run, created by this script

Each run folder holds:

    models/<name>/runs/<timestamp>/
        input/00001_input.txt     template values of every sample, all written before
                                  the first simulation runs
        output/00001_output.txt   PHREEQC selected output of that simulation
        failed/00001_input.txt    input files of failed, timed-out and flagged
                                  simulations, moved here from input/
        report.txt                configuration, seed, counts, ranges, failure reasons
        config.cfg                copy of the config used
        <name>.phr                copy of the template used
        run.log, logs/, plot/     only when debug = true
"""

import configparser
import hashlib
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
from datetime import datetime
from multiprocessing import Pool
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import qmc, triang, truncnorm

log = logging.getLogger("datagen")


class SimulationError(Exception):
    """A run finished but its result cannot be used."""


# =============================================================================
# Setup
# =============================================================================


def items(value):
    """Comma-separated list from the config; empty value = empty list."""
    return [x.strip() for x in value.split(",") if x.strip()]


def load_config(path):
    """Read config.cfg into the nested structure the rest of the script uses."""
    cp = configparser.ConfigParser(interpolation=None, inline_comment_prefixes=("#",))
    cp.optionxform = str  # Keep key case as written.
    with open(path, encoding="utf-8") as f:
        cp.read_file(f)

    run, ph, smp, model, chk, out = (
        cp[s] for s in ("run", "phreeqc", "sampling", "model", "checks", "output")
    )
    text_keys = {
        "distribution",
        "placeholder",
        "to",
        "fraction_parameter",
        "moles_placeholder",
        "area_placeholder",
    }

    def numbers(section):
        return {k: v if k in text_keys else float(v) for k, v in section.items()}

    return {
        "run": {
            "enabled": run.getboolean("enabled"),
            "model_encoding": run["model_encoding"],
            "n_simulations": run.getint("n_simulations"),
            "debug": run.getboolean("debug"),
            "workers": run.getint("workers"),
            "timeout_s": run.getfloat("timeout_s"),
            "smoke_test": run.getboolean("smoke_test"),
        },
        "phreeqc": {"executable": ph["executable"], "database": ph["database"]},
        "sampling": {
            "method": smp["method"],
            "seed": smp["seed"] if smp["seed"] == "auto" else int(smp["seed"]),
        },
        "model": {
            "volume_m3": model.getfloat("volume_m3"),
            "total_pressure_atm": model.getfloat("total_pressure_atm"),
            "h2_pressure_atm": model.getfloat("h2_pressure_atm"),
            "normalize_mineral_fractions": model.getboolean(
                "normalize_mineral_fractions"
            ),
        },
        "parameters": {
            name.split(".", 1)[1]: numbers(cp[name])
            for name in cp.sections()
            if name.startswith("parameter.")
        },
        "minerals": [
            {"name": name.split(".", 1)[1], **numbers(cp[name])}
            for name in cp.sections()
            if name.startswith("mineral.")
        ],
        "checks": {
            "time_column": chk["time_column"],
            "end_time_tolerance": chk.getfloat("end_time_tolerance"),
            "ph_column": chk["ph_column"],
            "ph_range": [float(x) for x in items(chk["ph_range"])],
            "nonnegative": items(chk["nonnegative"]),
            "must_change": items(chk["must_change"]),
            "change_atol": chk.getfloat("change_atol"),
        },
        "output": {
            "input_file": items(out["input_file"]),
            "drop_columns": items(out["drop_columns"]),
            "include_flagged": out.getboolean("include_flagged"),
            "keep_initial_row": out.getboolean("keep_initial_row"),
        },
        "plot": [
            {
                "title": name.split(".", 1)[1],
                "columns": items(cp[name]["columns"]),
                "scale": cp[name].getfloat("scale"),
                "ylabel": cp[name]["ylabel"],
            }
            for name in cp.sections()
            if name.startswith("plot.")
        ],
    }


def setup_logging(log_file):
    log.handlers.clear()  # One model after another: start each with fresh handlers.
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.INFO)
    console.setFormatter(
        logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%H:%M:%S")
    )
    log.addHandler(console)
    if log_file is not None:
        fh = logging.FileHandler(log_file, encoding="utf-8")
        fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)-7s %(message)s"))
        log.addHandler(fh)


# =============================================================================
# Sampling
# =============================================================================


def unit_points(method, d, seed, n):
    """n points in [0, 1)^d."""
    if method == "lhs":
        return qmc.LatinHypercube(d=d, seed=seed).random(n)
    if method == "sobol":
        return qmc.Sobol(d=d, scramble=True, seed=seed).random(n)
    return np.random.default_rng(seed).random((n, d))


# Each distribution maps a unit-cube value u to a parameter value through its
# inverse CDF, so Latin hypercube stratification carries over. All of them are
# bounded by low and high (normal and lognormal are truncated there).


def _truncated_normal(u, mean, std, lo, hi):
    return truncnorm.ppf(u, (lo - mean) / std, (hi - mean) / std, loc=mean, scale=std)


def _uniform(u, p):
    return p["low"] + u * (p["high"] - p["low"])


def _loguniform(u, p):
    lo, hi = np.log10(p["low"]), np.log10(p["high"])
    return 10 ** (lo + u * (hi - lo))


def _triangular(u, p):
    width = p["high"] - p["low"]
    return triang.ppf(u, (p["mode"] - p["low"]) / width, loc=p["low"], scale=width)


def _normal(u, p):
    return _truncated_normal(u, p["mean"], p["std"], p["low"], p["high"])


def _lognormal(u, p):
    return 10 ** _truncated_normal(
        u,
        np.log10(p["median"]),
        p["sigma_log10"],
        np.log10(p["low"]),
        np.log10(p["high"]),
    )


DISTRIBUTIONS = {
    "uniform": _uniform,  # low, high
    "loguniform": _loguniform,  # low, high
    "triangular": _triangular,  # low, mode, high
    "normal": _normal,  # mean, std; truncated to [low, high]
    "lognormal": _lognormal,  # median, sigma_log10; truncated to [low, high]
}

# Parameters that are not sampled and use no dimension of the unit cube:
#   constant      value
#   proportional  to, base, base_to: value = sample[to] * base / base_to
COMPUTED = ("constant", "proportional")


def sampled(cfg):
    """Names of the parameters drawn from a distribution, in config order."""
    return [
        n for n, p in cfg["parameters"].items() if p["distribution"] not in COMPUTED
    ]


def scale(unit_row, cfg):
    params = cfg["parameters"]
    s = {
        name: float(DISTRIBUTIONS[params[name]["distribution"]](u, params[name]))
        for u, name in zip(unit_row, sampled(cfg))
    }
    for name, p in params.items():
        if p["distribution"] == "constant":
            s[name] = p["value"]
    for name, p in params.items():
        if p["distribution"] == "proportional":
            s[name] = s[p["to"]] * p["base"] / p["base_to"]
    return {name: s[name] for name in params}


def derive(s, cfg):
    """Quantities computed from the sampled parameters."""
    m = cfg["model"]
    d = {}

    pore_volume = m["volume_m3"] * s["porosity"] * 1000.0  # L
    d["water_volume_L"] = pore_volume * s["water_fraction"]
    d["gas_volume_L"] = pore_volume - d["water_volume_L"]

    d["pCH4"] = (
        m["total_pressure_atm"]
        - s["pCO2"]
        - s["pH2S"]
        - s["pN2"]
        - m["h2_pressure_atm"]
    )

    # Volume fraction of the solid for each mineral, optionally scaled to sum to 1.
    fractions = {mn["name"]: s[mn["fraction_parameter"]] for mn in cfg["minerals"]}
    if m["normalize_mineral_fractions"]:
        total = sum(fractions.values())
        fractions = {name: f / total for name, f in fractions.items()}

    for mineral in cfg["minerals"]:
        name = mineral["name"]
        mass = (
            mineral["density_kg_m3"]
            * m["volume_m3"]
            * (1 - s["porosity"])
            * fractions[name]
        )
        rf = mineral["reactive_fraction"]
        d[f"{name}_area_m2"] = (
            (6 / mineral["diameter_m"]) * (mass / mineral["density_kg_m3"]) * rf
        ) * mineral["area_multiplier"]
        d[f"{name}_moles"] = (
            mass / mineral["molar_mass_kg_mol"] * rf * mineral["moles_multiplier"]
        )
    return d


def replacements(s, d, cfg):
    """Placeholder -> value text, formatted as in the original script (.4e).

    The same text goes into the template and into the input file, so the
    input file holds exactly what PHREEQC was given.
    """
    reps = {
        p["placeholder"]: f"{s[name]:.4e}"
        for name, p in cfg["parameters"].items()
        if "placeholder" in p
    }
    reps.update(
        {
            "{WATER_VOLUME}": f"{d['water_volume_L']:.4e}",
            "{GAS_VOLUME}": f"{d['gas_volume_L']:.4e}",
            "{CH4}": f"{d['pCH4']:.4e}",
        }
    )
    for mineral in cfg["minerals"]:
        name = mineral["name"]
        reps[mineral["moles_placeholder"]] = f"{d[f'{name}_moles']:.4e}"
        reps[mineral["area_placeholder"]] = f"{d[f'{name}_area_m2']:.4e}"
    return reps


def render(template, reps):
    for key, value in reps.items():
        template = template.replace(key, value)
    return template


def end_days(template):
    """Total kinetic time from the '-steps <seconds> in <n>' line, in days."""
    return (
        float(re.search(r"^\s*-steps\s+([0-9.eE+]+)", template, re.MULTILINE).group(1))
        / 86400.0
    )


# =============================================================================
# Worker side
# =============================================================================

_W = {}  # Per-process settings, set once by init_worker.


def init_worker(settings):
    _W.update(settings)


def run_phreeqc(text, tag):
    logs = Path(_W["run_dir"]) / "logs"
    with tempfile.TemporaryDirectory() as tmp:
        pqi, pqo = os.path.join(tmp, "model.pqi"), os.path.join(tmp, "model.pqo")
        Path(pqi).write_text(text, encoding=_W["encoding"])
        if _W["debug"]:
            shutil.copyfile(pqi, logs / f"{tag}_input.pqi")

        try:
            res = subprocess.run(
                [_W["phreeqc"], pqi, pqo, _W["database"]],
                capture_output=True,
                text=True,
                cwd=tmp,
                timeout=_W["timeout"],
                check=False,
            )
        except subprocess.TimeoutExpired:
            raise SimulationError(f"timeout after {_W['timeout']} s")

        out = Path(pqo).read_text(encoding=_W["encoding"])
        if _W["debug"]:
            (logs / f"{tag}_output.pqo").write_text(
                out + res.stderr, encoding=_W["encoding"]
            )

        errors = [
            ln.strip()
            for ln in out.splitlines()
            if ln.startswith("ERROR") and ln.strip() != "ERROR:"
        ]
        if errors:
            raise SimulationError(" | ".join(errors[:3]))
        if res.returncode != 0:
            raise SimulationError(f"PHREEQC exit code {res.returncode}")
        return Path(tmp, "selected_output_1.sel").read_text(encoding=_W["encoding"])


def first_step(time_column):
    """Index of the first kinetic step (t > 0).

    Rows before it come from PHREEQC's initial solution calculation, done
    before the gas phase and kinetics are set up (default pH 7, no gas,
    KIN() = 0). They are not a physical state of the system.
    """
    return int(np.argmax(time_column > 0))


def quality_issues(col, data):
    """Problems that make a completed run suspicious. Empty list = clean."""
    chk = _W["checks"]
    issues = []

    ph = data[:, col[chk["ph_column"]]]
    lo, hi = chk["ph_range"]
    if ph.min() < lo or ph.max() > hi:
        issues.append(f"pH outside [{lo}, {hi}] ({ph.min():.2f} to {ph.max():.2f})")

    negative = [c for c in chk["nonnegative"] if data[:, col[c]].min() < 0]
    if negative:
        issues.append(f"negative values in {', '.join(negative)}")

    # Absolute change: kinetic reactants often start at an arbitrary -m 1,
    # so real progress can be tiny relative to the value.
    unchanged = [
        c
        for c in chk["must_change"]
        if abs(data[-1, col[c]] - data[0, col[c]]) <= chk["change_atol"]
    ]
    if unchanged:
        issues.append(f"no change in {', '.join(unchanged)}")
    return issues


def work(task):
    """Run one simulation. Returns a result dict; errors become a failed status."""
    sim_id, tag = task["sim_id"], task["tag"]
    t0 = time.perf_counter()
    res = {
        "sim_id": sim_id,
        "status": "failed",
        "message": "",
        "headers": None,
        "data": None,
    }
    try:
        lines = run_phreeqc(render(_W["template"], task["reps"]), tag).splitlines()
        headers = lines[0].split()
        data = np.array([[float(x) for x in ln.split()] for ln in lines[1:]])
        col = {h: i for i, h in enumerate(headers)}
        chk = _W["checks"]

        t = data[:, col[chk["time_column"]]]
        if t[-1] < chk["end_time_tolerance"] * _W["end_days"]:
            raise SimulationError(
                f"stopped early at {t[-1]:.2f} of {_W['end_days']:.2f} d"
            )

        kinetic = data[first_step(t) :]
        issues = quality_issues(col, kinetic)
        res["status"] = "flagged" if issues else "success"
        res["message"] = "; ".join(issues)
        res["headers"] = headers
        res["data"] = data if _W["keep_initial_row"] else kinetic

        if _W["debug"]:
            plot_run(
                sim_id, col, kinetic, Path(_W["run_dir"]) / "plot" / f"{tag}_plot.jpg"
            )

    except SimulationError as e:
        res["status"] = "timeout" if str(e).startswith("timeout") else "failed"
        res["message"] = str(e)
    except Exception as e:
        res["status"], res["message"] = "failed", f"{type(e).__name__}: {e}"
        if _W["debug"]:
            (Path(_W["run_dir"]) / "logs" / f"{tag}_traceback.txt").write_text(
                traceback.format_exc(), encoding="utf-8"
            )

    res["runtime_s"] = round(time.perf_counter() - t0, 2)
    return res


def plot_run(sim_id, col, data, path):
    panels = _W["plot"]
    t = data[:, col[_W["checks"]["time_column"]]]

    nrows = (len(panels) + 1) // 2
    fig, axes = plt.subplots(
        nrows=nrows, ncols=2, figsize=(14, 5 * nrows), squeeze=False
    )
    fig.suptitle(f"Simulation {sim_id}", fontsize=14, fontweight="bold")
    for ax, panel in zip(axes.flat, panels):
        for c in panel["columns"]:
            ax.plot(t, data[:, col[c]] * panel["scale"], label=c)
        ax.set_title(panel["title"])
        ax.set_xlabel("Time (days)")
        ax.set_ylabel(panel["ylabel"])
        ax.grid(True, linestyle="--", alpha=0.7)
        if len(panel["columns"]) > 1:
            ax.legend(loc="best", fontsize="small", ncol=2)
    for ax in list(axes.flat)[len(panels) :]:
        ax.set_visible(False)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=150)
    plt.close(fig)


# =============================================================================
# Result files
# =============================================================================


def write_input(path, reps, keys):
    """Same layout as the original script: one '{PLACEHOLDER} value' line each."""
    path.write_text("".join(f"{k} {reps[k]}\n" for k in keys), encoding="utf-8")


def write_output(path, headers, data, drop):
    keep = [i for i, h in enumerate(headers) if h not in drop]
    lines = ["\t".join(headers[i] for i in keep)]
    lines += ["\t".join(f"{v:.12e}" for v in line[keep]) for line in data]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def save_result(run_dir, task, res, cfg, usable):
    """Write the output of a usable run; move the input file of any other run."""
    tag = task["tag"]
    if res["status"] in usable:
        write_output(
            run_dir / "output" / f"{tag}_output.txt",
            res["headers"],
            res["data"],
            set(cfg["output"]["drop_columns"]),
        )
    else:
        shutil.move(
            run_dir / "input" / f"{tag}_input.txt",
            run_dir / "failed" / f"{tag}_input.txt",
        )


# =============================================================================
# Report
# =============================================================================


def write_report(run_dir, cfg, meta, records, n_planned):
    counts = {s: 0 for s in ("success", "flagged", "failed", "timeout")}
    for _, res in records:
        counts[res["status"]] += 1
    sep = "=" * 80
    out = [
        sep,
        "General Information",
        sep,
        "",
        f"Run folder      => {run_dir}",
        f"Model           => {meta['template_name']} (sha256 {meta['template_sha256'][:12]})",
        f"PHREEQC         => {cfg['phreeqc']['executable']}",
        f"Database        => {cfg['phreeqc']['database']}",
        f"Started         => {meta['started']}",
        f"Mode            => {'debug' if cfg['run']['debug'] else 'production'}",
        f"Sampling        => {cfg['sampling']['method']}, seed {meta['seed']}",
        f"Model encoding  => {cfg['run']['model_encoding']}",
        f"Workers         => {cfg['run']['workers']}",
        f"Timeout         => {cfg['run']['timeout_s']} s",
        f"Smoke test      => {meta['smoke']}",
        f"Simulations     => {len(records)} of {n_planned} finished",
    ]
    out += [f"  {status:<14}=> {count}" for status, count in counts.items()]
    if records:
        runtimes = [res["runtime_s"] for _, res in records]
        out.append(
            f"Runtime per sim => median {np.median(runtimes):.1f} s, max {max(runtimes):.1f} s"
        )

    ok = [row for row, res in records if res["status"] == "success"]
    bad = [row for row, res in records if res["status"] != "success"]

    def rng(rows, name):
        v = [r[name] for r in rows]
        return f"{min(v):>12.4e}{max(v):>12.4e}" if v else f"{'-':>12}{'-':>12}"

    out += [
        "",
        sep,
        "Parameter Ranges",
        sep,
        "",
        f"  {'parameter':<16}{'distribution':<12}{'bounds':>24}  {'successful runs':>24}  {'other runs':>24}",
    ]
    for name, p in cfg["parameters"].items():
        bounds = (
            f"{p['low']:>12.4e}{p['high']:>12.4e}"
            if p["distribution"] not in COMPUTED
            else f"{'-':>12}{'-':>12}"
        )
        out.append(
            f"  {name:<16}{p['distribution']:<12}{bounds}  "
            f"{rng(ok, name)}  {rng(bad, name)}"
        )

    problems = sorted(
        [(row, res) for row, res in records if res["status"] != "success"],
        key=lambda r: r[0]["sim_id"],
    )
    if problems:
        out += ["", sep, "Flagged and Failed Runs", sep, ""]
        out += [
            f"  {row['sim_id']:>6}  {res['status']:<8} {res['message']}"
            for row, res in problems
        ]
    (run_dir / "report.txt").write_text("\n".join(out) + "\n", encoding="utf-8")


# =============================================================================
# Main
# =============================================================================


MODELS_DIR = Path(__file__).resolve().parent / "models"


def run_model(model_dir, cfg):
    """Sample, run and save one model. Returns False if its smoke test fails."""
    run = cfg["run"]
    model_path = model_dir / f"{model_dir.name}.phr"
    template = model_path.read_text(encoding=run["model_encoding"])

    seed = cfg["sampling"]["seed"]
    if seed == "auto":
        seed = int(np.random.SeedSequence().entropy % 2**32)
    n = run["n_simulations"]

    # All samples are drawn up front in this process from one seed.
    rows = []
    for i, u in enumerate(
        unit_points(cfg["sampling"]["method"], len(sampled(cfg)), seed, n), 1
    ):
        s = scale(u, cfg)
        rows.append({"sim_id": i, **s, **derive(s, cfg)})

    stamp = datetime.now()
    run_dir = model_dir / "runs" / f"{stamp:%Y%m%d-%H%M%S}"
    folders = ["input", "output", "failed"] + (["logs", "plot"] if run["debug"] else [])
    for folder in folders:
        (run_dir / folder).mkdir(parents=True)
    shutil.copyfile(model_dir / "config.cfg", run_dir / "config.cfg")
    shutil.copyfile(model_path, run_dir / model_path.name)
    setup_logging(run_dir / "run.log" if run["debug"] else None)

    log.info("Model %s | run folder: %s", model_dir.name, run_dir)
    log.info(
        "Mode: %s | seed %d | %d samples (%s)",
        "debug" if run["debug"] else "production",
        seed,
        n,
        cfg["sampling"]["method"],
    )

    meta = {
        "seed": seed,
        "started": stamp.isoformat(timespec="seconds"),
        "template_name": model_path.name,
        "template_sha256": hashlib.sha256(model_path.read_bytes()).hexdigest(),
        "smoke": "on" if run["smoke_test"] else "off",
    }
    settings = {
        "template": template,
        "phreeqc": cfg["phreeqc"]["executable"],
        "database": cfg["phreeqc"]["database"],
        "timeout": run["timeout_s"],
        "debug": run["debug"],
        "run_dir": str(run_dir),
        "checks": cfg["checks"],
        "end_days": end_days(template),
        "plot": cfg["plot"],
        "keep_initial_row": cfg["output"]["keep_initial_row"],
        "encoding": run["model_encoding"],
    }
    usable = (
        ("success", "flagged") if cfg["output"]["include_flagged"] else ("success",)
    )

    # Every input file is written before anything runs, so the whole sample set
    # from this seed is on disk from the start.
    tasks = [
        {
            "sim_id": r["sim_id"],
            "tag": f"{r['sim_id']:05d}",
            "reps": replacements({k: r[k] for k in cfg["parameters"]}, r, cfg),
        }
        for r in rows
    ]
    for t in tasks:
        write_input(
            run_dir / "input" / f"{t['tag']}_input.txt",
            t["reps"],
            cfg["output"]["input_file"],
        )
    log.info("Wrote %d input files to %s", n, run_dir / "input")
    write_report(run_dir, cfg, meta, [], n)  # Seed and settings on disk from the start.

    if run["smoke_test"]:
        # The median of every parameter distribution, run once before the batch.
        s = scale([0.5] * len(sampled(cfg)), cfg)
        init_worker(settings)
        log.info("Smoke test: running the median of every parameter...")
        res = work(
            {"sim_id": 0, "tag": "00000", "reps": replacements(s, derive(s, cfg), cfg)}
        )
        meta["smoke"] = (
            f"{res['status']} in {res['runtime_s']} s {res['message']}".strip()
        )
        if res["status"] != "success":
            log.error("Smoke test %s: %s", res["status"], res["message"])
            log.error("Skipping %s.", model_dir.name)
            write_report(run_dir, cfg, meta, [], n)
            return False
        log.info(
            "Smoke test passed in %.1f s. Rough batch estimate: %.1f min.",
            res["runtime_s"],
            res["runtime_s"] * n / run["workers"] / 60,
        )

    by_id = {r["sim_id"]: r for r in rows}
    task_by_id = {t["sim_id"]: t for t in tasks}
    log.info(
        "Running %d simulations on %d workers (timeout: %s s).",
        n,
        run["workers"],
        run["timeout_s"],
    )

    t0 = time.perf_counter()
    records = []
    with Pool(
        processes=run["workers"], initializer=init_worker, initargs=(settings,)
    ) as pool:
        for res in pool.imap_unordered(work, tasks):
            row = by_id[res["sim_id"]]
            save_result(run_dir, task_by_id[res["sim_id"]], res, cfg, usable)
            records.append((row, res))
            level = logging.INFO if res["status"] == "success" else logging.WARNING
            note = f" ({res['message']})" if res["message"] else ""
            log.log(
                level,
                "[%d/%d] sim %d: %s in %.1f s%s",
                len(records),
                n,
                res["sim_id"],
                res["status"],
                res["runtime_s"],
                note,
            )

    write_report(run_dir, cfg, meta, records, n)
    n_ok = sum(res["status"] in usable for _, res in records)
    log.info(
        "Finished in %.1f s: %d of %d simulations in output/.",
        time.perf_counter() - t0,
        n_ok,
        n,
    )
    log.info("Report: %s", run_dir / "report.txt")
    return True


def main():
    model_dirs = sorted(d for d in MODELS_DIR.iterdir() if (d / "config.cfg").is_file())
    results = {}
    for model_dir in model_dirs:
        cfg = load_config(model_dir / "config.cfg")
        if cfg["run"]["enabled"]:
            results[model_dir.name] = run_model(model_dir, cfg)
        else:
            results[model_dir.name] = None

    setup_logging(None)
    for name, ok in results.items():
        status = {True: "done", False: "smoke test failed", None: "disabled"}[ok]
        log.info("%-30s %s", name, status)


if __name__ == "__main__":
    main()
