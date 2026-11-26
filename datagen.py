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
from datetime import datetime, timedelta, timezone
from multiprocessing import Pool
from pathlib import Path
from statistics import NormalDist

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

log = logging.getLogger("datagen")


class SimulationError(Exception):
    """A run finished but its result cannot be used."""


def items(value):
    """Comma-separated list from the config; empty value = empty list."""
    return [x.strip() for x in value.split(",") if x.strip()]


def load_config(path):
    """Read config.cfg into the nested structure the rest of the script uses."""
    cp = configparser.ConfigParser(interpolation=None, inline_comment_prefixes=("#",))
    cp.optionxform = str
    with open(path, encoding="utf-8") as f:
        cp.read_file(f)

    run, ph, model, chk, out = (
        cp[s] for s in ("run", "phreeqc", "model", "checks", "output")
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

    if model.get("cushion_gas") not in GASES:
        raise ValueError(f"{path}: cushion_gas must be one of {', '.join(GASES)}")
    known = [*DISTRIBUTIONS, *COMPUTED]
    for name in cp.sections():
        if name.startswith("parameter.") and cp[name]["distribution"] not in known:
            raise ValueError(
                f"{path}: [{name}] distribution must be one of {', '.join(known)}"
            )

    return {
        "run": {
            "enabled": run.getboolean("enabled"),
            "n_simulations": run.getint("n_simulations"),
            "debug": run.getboolean("debug"),
            "workers": run.getint("workers"),
            "timeout_s": run.getfloat("timeout_s"),
        },
        "phreeqc": {"executable": ph["executable"], "database": ph["database"]},
        "model": {
            "cushion_gas": model["cushion_gas"],
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
    log.handlers.clear()
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


def latin_hypercube(d, seed, n):
    """n points in [0, 1)^d, one in each of n equal slices of every dimension.

    The same algorithm and random calls as scipy.stats.qmc.LatinHypercube
    (scramble=True), so a seed gives the same points as with scipy.
    """
    rng = np.random.default_rng(seed)
    offsets = rng.uniform(size=(n, d))
    slices = np.tile(np.arange(1, n + 1), (d, 1))
    for i in range(d):
        rng.shuffle(slices[i, :])
    return (slices.T - offsets) / n


_STANDARD_NORMAL = NormalDist()


def _truncated_normal(u, mean, std, lo, hi):
    """Inverse CDF of a normal distribution truncated to [lo, hi]."""
    a = _STANDARD_NORMAL.cdf((lo - mean) / std)
    b = _STANDARD_NORMAL.cdf((hi - mean) / std)
    return mean + std * _STANDARD_NORMAL.inv_cdf(a + u * (b - a))


def _uniform(u, p):
    return p["low"] + u * (p["high"] - p["low"])


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
    "uniform": _uniform,
    "normal": _normal,
    "lognormal": _lognormal,
}

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


GASES = ("CO2", "H2S", "N2", "CH4")

# Every model is a 1 m³ REV of rock.
VOLUME_M3 = 1.0

# Templates and every file datagen.py writes are UTF-8.
ENCODING = "utf-8"


def gas_pressures(template):
    """Total and H2 pressure (atm) from the template's GAS_PHASE.

    Read from the "-pressure" and "Hdg(g)" lines, so the cushion gas always
    takes exactly what the template leaves of the total.
    """

    def value(pattern, what):
        found = re.search(pattern, template, re.MULTILINE)
        if not found:
            raise ValueError(f"template has no {what} line in GAS_PHASE")
        return float(found.group(1))

    number = r"([0-9.]+(?:[eE][+-]?[0-9]+)?)"
    return {
        "total_pressure_atm": value(rf"^\s*-pressure\s+{number}", '"-pressure"'),
        "h2_pressure_atm": value(rf"^\s*Hdg\(g\)\s+{number}", '"Hdg(g)"'),
    }


def derive(s, cfg):
    """Quantities computed from the sampled parameters."""
    m = cfg["model"]
    d = {}

    pore_volume = VOLUME_M3 * s["porosity"] * 1000.0
    d["water_volume_L"] = pore_volume * s["water_fraction"]
    d["gas_volume_L"] = pore_volume - d["water_volume_L"]

    cushion = m["cushion_gas"]
    d[f"p{cushion}"] = (
        m["total_pressure_atm"]
        - m["h2_pressure_atm"]
        - sum(s[f"p{gas}"] for gas in GASES if gas != cushion)
    )

    # The mineral parameters are mass shares of the solid (e.g. percent).
    # Volume share = mass share / density, scaled to sum to 1.
    shares = {mn["name"]: s[mn["fraction_parameter"]] for mn in cfg["minerals"]}
    total = sum(shares.values())
    d.update({f"{n}_mass_fraction": f / total for n, f in shares.items()})
    volume = {
        mn["name"]: shares[mn["name"]] / mn["density_kg_m3"] for mn in cfg["minerals"]
    }
    total = sum(volume.values())
    fractions = {name: v / total for name, v in volume.items()}

    for mineral in cfg["minerals"]:
        name = mineral["name"]
        d[f"{name}_volume_fraction"] = fractions[name]
        mass = (
            mineral["density_kg_m3"] * VOLUME_M3 * (1 - s["porosity"]) * fractions[name]
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
        }
    )
    cushion = cfg["model"]["cushion_gas"]
    reps[f"{{{cushion}}}"] = f"{d[f'p{cushion}']:.4e}"
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


_W = {}


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
        version = re.search(
            r"\*\s*(PHREEQC[\w.\-]*)\s*\*", res.stdout + res.stderr + out
        )
        return (
            Path(tmp, "selected_output_1.sel").read_text(encoding=_W["encoding"]),
            version.group(1) if version else None,
        )


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
        "phreeqc_version": None,
    }
    try:
        sel, res["phreeqc_version"] = run_phreeqc(
            render(_W["template"], task["reps"]), tag
        )
        lines = sel.splitlines()
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
    except (OSError, ValueError, KeyError, IndexError) as e:
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
        res["n_rows"] = len(res["data"])
    else:
        shutil.move(
            run_dir / "input" / f"{tag}_input.txt",
            run_dir / "failed" / f"{tag}_input.txt",
        )


STATUSES = ("success", "flagged", "failed", "timeout")
SEP = "=" * 80


def _num(v):
    return f"{v:.4e}"


def _table(header, lines):
    """Text table: first column left-aligned, the rest right-aligned."""
    widths = [max(len(str(c)) for c in col) for col in zip(header, *lines)]

    def fmt(cells):
        return "  " + "  ".join(
            f"{c:<{w}}" if k == 0 else f"{c:>{w}}"
            for k, (c, w) in enumerate(zip(map(str, cells), widths))
        )

    return [fmt(header), "  " + "  ".join("-" * w for w in widths)] + [
        fmt(line) for line in lines
    ]


def _stats(rows, name):
    """min, median, max of one quantity over some runs, as text."""
    v = [r[name] for r in rows]
    if not v:
        return ["-", "-", "-"]
    return [_num(min(v)), _num(float(np.median(v))), _num(max(v))]


def _id_ranges(ids):
    """[1, 2, 3, 7] -> '1-3, 7'."""
    ids, parts = sorted(ids), []
    start = prev = ids[0]
    for k in ids[1:] + [None]:
        if k is not None and k == prev + 1:
            prev = k
            continue
        parts.append(f"{start}" if start == prev else f"{start}-{prev}")
        if k is not None:
            start = prev = k
    return ", ".join(parts)


def _section(title):
    return ["", SEP, title, SEP, ""]


def write_report(run_dir, cfg, meta, records, n_planned):
    """Write report.txt: settings, results, value ranges, failures and every run.

    Called before the batch, during it and at the end, so the report is
    always current, even if the batch is stopped.
    """
    params = cfg["parameters"]
    varied = [n for n, p in params.items() if p["distribution"] != "constant"]
    derived = meta["derived_keys"]
    fractions = [
        k for k in derived if k.endswith(("_mass_fraction", "_volume_fraction"))
    ]
    counts = {s: sum(res["status"] == s for _, res in records) for s in STATUSES}
    ok = [row for row, res in records if res["status"] == "success"]
    other = [row for row, res in records if res["status"] != "success"]
    m = cfg["model"]

    out = [SEP, "General Information", SEP, ""]
    out += [
        f"Run folder      => {run_dir}",
        f"Model           => {meta['template_name']} (sha256 {meta['template_sha256'][:12]})",
        f"datagen.py      => {meta['code_version']}",
        f"PHREEQC         => {cfg['phreeqc']['executable']} ({meta.get('phreeqc_version') or 'version not known yet'})",
        f"Database        => {cfg['phreeqc']['database']}",
        f"Started         => {meta['started']}",
        f"Finished        => {meta['finished'] or 'still running'}",
        f"Mode            => {'debug' if cfg['run']['debug'] else 'production'}",
        f"Sampling        => Latin hypercube, seed {meta['seed']}",
        f"Workers         => {cfg['run']['workers']}",
        f"Timeout         => {cfg['run']['timeout_s']} s",
        f"Gas             => {m['total_pressure_atm']} atm total, "
        f"H2 {m['h2_pressure_atm']} atm, cushion gas {m['cushion_gas']}",
        f"Simulations     => {len(records)} of {n_planned} finished",
    ]
    out += [f"  {status:<14}=> {count}" for status, count in counts.items()]
    if records:
        runtimes = [res["runtime_s"] for _, res in records]
        out.append(
            f"Runtime per sim => median {np.median(runtimes):.1f} s, "
            f"max {max(runtimes):.1f} s, total {sum(runtimes) / 3600:.2f} CPU h"
        )
    saved = [res for _, res in records if res.get("n_rows")]
    if saved:
        n_rows = [res["n_rows"] for res in saved]
        rows_text = (
            f"{n_rows[0]} rows each"
            if min(n_rows) == max(n_rows)
            else f"{min(n_rows)} to {max(n_rows)} rows"
        )
        out.append(
            f"Output          => {len(saved)} files in output/, {rows_text}, "
            f"{len(meta['output_columns'])} columns"
        )

    out += _section("Sampled Parameters")
    out += [
        "Ranges over the runs, split into successful runs and all others",
        "(flagged, failed, timeout). Constant parameters are listed once.",
        "",
    ]
    lines = []
    for name, p in params.items():
        if p["distribution"] == "constant":
            lines.append(
                [name, "constant", "-", "-", _num(p["value"]), "", "", "", "", ""]
            )
            continue
        bounds = (
            [_num(p["low"]), _num(p["high"])]
            if p["distribution"] not in COMPUTED
            else ["-", "-"]
        )
        lines.append(
            [name, p["distribution"], *bounds, *_stats(ok, name), *_stats(other, name)]
        )
    out += _table(
        [
            "parameter",
            "distribution",
            "low",
            "high",
            "ok min",
            "ok median",
            "ok max",
            "other min",
            "other median",
            "other max",
        ],
        lines,
    )

    out += _section("Derived Values")
    out += [
        "Computed from the parameters (see the config): volumes in L, gas",
        "pressure in atm, mineral moles and areas in m², and each mineral's",
        "mass and volume fractions of the solid after scaling to sum to 1.",
        "",
    ]
    out += _table(
        [
            "quantity",
            "ok min",
            "ok median",
            "ok max",
            "other min",
            "other median",
            "other max",
        ],
        [[k, *_stats(ok, k), *_stats(other, k)] for k in derived],
    )

    problems = [(row, res) for row, res in records if res["status"] != "success"]
    if problems:
        out += _section("Flagged and Failed Runs")
        groups = {}
        for row, res in problems:
            groups.setdefault((res["status"], res["message"]), []).append(row["sim_id"])
        for (status, message), ids in sorted(groups.items(), key=lambda g: -len(g[1])):
            out += [
                f"{status}: {len(ids)} run{'s' if len(ids) > 1 else ''}",
                f"  {message}",
                f"  runs: {_id_ranges(ids)}",
                "",
            ]

    if records:
        out += _section("Simulations")
        out += [
            "One line per finished run: the varied parameters, with the mineral",
            "fractions as used (after scaling to sum to 1). The values PHREEQC",
            "was given are in input/ (or failed/) under the same number.",
            "",
        ]
        drawn = {mn["fraction_parameter"] for mn in cfg["minerals"]}
        cols = [n for n in varied if n not in drawn] + fractions
        out += _table(
            ["sim", "status", "runtime_s", *cols],
            [
                [
                    row["sim_id"],
                    res["status"],
                    f"{res['runtime_s']:.1f}",
                    *[_num(row[c]) for c in cols],
                ]
                for row, res in sorted(records, key=lambda r: r[0]["sim_id"])
            ],
        )

    if meta.get("output_columns"):
        out += _section("Output Columns")
        out += [", ".join(meta["output_columns"])]

    (run_dir / "report.txt").write_text("\n".join(out) + "\n", encoding="utf-8")


def code_version():
    """sha256 of datagen.py, plus its git commit when this is a git checkout."""
    path = Path(__file__).resolve()
    text = f"sha256 {hashlib.sha256(path.read_bytes()).hexdigest()[:12]}"
    try:
        git = ["git", "-C", str(path.parent)]
        commit = subprocess.run(
            git + ["rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
            timeout=10,
        ).stdout.strip()
        dirty = subprocess.run(
            git + ["status", "--porcelain", "--", path.name],
            capture_output=True,
            text=True,
            check=True,
            timeout=10,
        ).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return text
    return f"{text}, git {commit}" + (" with uncommitted changes" if dirty else "")


MODELS_DIR = Path(__file__).resolve().parent / "models"

# Run folders are named in Istanbul time. Turkey has kept UTC+3 all year
# since 2016, so a fixed offset needs no time-zone database.
ISTANBUL = timezone(timedelta(hours=3), "Europe/Istanbul")


def now():
    return datetime.now(ISTANBUL)


def draw_samples(cfg, seed):
    """One row per simulation: its parameters and the values derived from them."""
    rows = []
    points = latin_hypercube(len(sampled(cfg)), seed, cfg["run"]["n_simulations"])
    for i, u in enumerate(points, 1):
        s = scale(u, cfg)
        rows.append({"sim_id": i, **s, **derive(s, cfg)})
    return rows


def make_run_dir(model_dir, model_path, debug, stamp):
    """Create the run folder with copies of the config and template."""
    run_dir = model_dir / "runs" / f"{stamp:%Y%m%d-%H%M%S}"
    for folder in ["input", "output", "failed"] + (["logs", "plot"] if debug else []):
        (run_dir / folder).mkdir(parents=True)
    shutil.copyfile(model_dir / "config.cfg", run_dir / "config.cfg")
    shutil.copyfile(model_path, run_dir / model_path.name)
    return run_dir


def write_inputs(run_dir, rows, cfg):
    """Write one input file per sample and return the tasks for the workers."""
    tasks = []
    for r in rows:
        task = {
            "sim_id": r["sim_id"],
            "tag": f"{r['sim_id']:05d}",
            "reps": replacements({k: r[k] for k in cfg["parameters"]}, r, cfg),
        }
        write_input(
            run_dir / "input" / f"{task['tag']}_input.txt",
            task["reps"],
            cfg["output"]["input_file"],
        )
        tasks.append(task)
    return tasks


def run_batch(run_dir, cfg, settings, meta, rows, tasks, usable):
    """Run every task on the worker pool, saving each result as it arrives."""
    n = len(tasks)
    by_id = {r["sim_id"]: r for r in rows}
    task_by_id = {t["sim_id"]: t for t in tasks}
    drop = set(cfg["output"]["drop_columns"])
    log.info(
        "Running %d simulations on %d workers (timeout: %s s).",
        n,
        cfg["run"]["workers"],
        cfg["run"]["timeout_s"],
    )

    last_report = time.perf_counter()
    records = []
    with Pool(
        processes=cfg["run"]["workers"], initializer=init_worker, initargs=(settings,)
    ) as pool:
        for res in pool.imap_unordered(work, tasks):
            save_result(run_dir, task_by_id[res["sim_id"]], res, cfg, usable)
            if res["headers"] and meta["output_columns"] is None:
                meta["output_columns"] = [h for h in res["headers"] if h not in drop]
            meta["phreeqc_version"] = meta["phreeqc_version"] or res["phreeqc_version"]
            res["data"] = None  # saved to output/; the report only needs the row count
            records.append((by_id[res["sim_id"]], res))

            note = f" ({res['message']})" if res["message"] else ""
            log.log(
                logging.INFO if res["status"] == "success" else logging.WARNING,
                "[%d/%d] sim %d: %s in %.1f s%s",
                len(records),
                n,
                res["sim_id"],
                res["status"],
                res["runtime_s"],
                note,
            )
            if time.perf_counter() - last_report > 60:
                write_report(run_dir, cfg, meta, records, n)
                last_report = time.perf_counter()
    return records


def run_model(model_dir, cfg):
    """Sample, run and save one model."""
    run = cfg["run"]
    model_path = model_dir / f"{model_dir.name}.phr"
    template = model_path.read_text(encoding=ENCODING)
    cfg["model"].update(gas_pressures(template))
    n = run["n_simulations"]

    seed = int(np.random.SeedSequence().entropy % 2**32)
    rows = draw_samples(cfg, seed)

    stamp = now()
    run_dir = make_run_dir(model_dir, model_path, run["debug"], stamp)
    setup_logging(run_dir / "run.log" if run["debug"] else None)
    log.info("Model %s | run folder: %s", model_dir.name, run_dir)
    log.info(
        "Mode: %s | seed %d | %d samples",
        "debug" if run["debug"] else "production",
        seed,
        n,
    )

    meta = {
        "seed": seed,
        "started": stamp.isoformat(timespec="seconds"),
        "finished": None,
        "template_name": model_path.name,
        "template_sha256": hashlib.sha256(model_path.read_bytes()).hexdigest(),
        "code_version": code_version(),
        "phreeqc_version": None,
        "output_columns": None,
        "derived_keys": [
            k for k in rows[0] if k != "sim_id" and k not in cfg["parameters"]
        ],
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
        "encoding": ENCODING,
    }
    usable = (
        ("success", "flagged") if cfg["output"]["include_flagged"] else ("success",)
    )

    tasks = write_inputs(run_dir, rows, cfg)
    log.info("Wrote %d input files to %s", n, run_dir / "input")
    write_report(run_dir, cfg, meta, [], n)

    t0 = time.perf_counter()
    records = run_batch(run_dir, cfg, settings, meta, rows, tasks, usable)
    wall = time.perf_counter() - t0
    meta["finished"] = (
        f"{now().isoformat(timespec='seconds')} (batch took {wall / 60:.1f} min)"
    )
    write_report(run_dir, cfg, meta, records, n)
    log.info(
        "Finished in %.1f s: %d of %d simulations in output/.",
        wall,
        sum(res["status"] in usable for _, res in records),
        n,
    )
    log.info("Report: %s", run_dir / "report.txt")


def main():
    model_dirs = sorted(d for d in MODELS_DIR.iterdir() if (d / "config.cfg").is_file())
    results = {}
    for model_dir in model_dirs:
        cfg = load_config(model_dir / "config.cfg")
        if cfg["run"]["enabled"]:
            run_model(model_dir, cfg)
            results[model_dir.name] = "done"
        else:
            results[model_dir.name] = "disabled"

    setup_logging(None)
    for name, status in results.items():
        log.info("%-30s %s", name, status)


if __name__ == "__main__":
    main()
