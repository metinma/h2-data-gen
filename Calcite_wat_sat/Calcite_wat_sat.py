import os
import shutil
import subprocess
import sys
import tempfile
from multiprocessing import Pool
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

MODEL = "Calcite_wat_sat.phr"
DEBUG = True
MAX_DURATION = 60


def main():
    WORKERS = max(1, (os.process_cpu_count() or 1) - 2)
    MODEL_SIZE = 100

    # Generate folders if they do not exist.
    generate_folders()

    args = [((i + 1), MODEL) for i in range(MODEL_SIZE)]

    # Run simulations with multiprocessing support.
    with Pool(processes=WORKERS) as pool:
        pool.starmap(work, args)

    print("All simulations have been completed.")

    # Carry the failed inputs another folder.
    organize_folders()

    # Print out the simulation report into a text file.
    generate_report(WORKERS, MODEL_SIZE)


def generate_folders():
    folders = ["input", "output"]
    if DEBUG:
        folders.append("plot")
        folders.append("logs")

    for folder in folders:
        os.makedirs(f"./{folder}", exist_ok=True)


# Normalizes a given list.
def normalize(bounds):
    factor = np.sum(bounds)
    normalized_bounds = bounds / factor
    return normalized_bounds


def original_values():
    solution_species = {
        "HCO3": 0.03445,
        "Na": 0.06278,
        "Mg": 0.00157,
        "Cl": 0.01162,
        "Ca": 0.00124,
        "S(6)": 0.01303,
    }

    gas_species = {
        "H2": 1.1843,
        "CH4": 53.4221,
        "CO2": 3.68428,
        "N2": 0.92107,
        "H2S": 0.003684,
    }

    # Fraction, Diameter (m), Density (kg/m^3), MW (kg/mol)
    mineral_species = {
        "Calcite": [1, 100e-6, 2710, 100.09e-3],
    }

    return solution_species, gas_species, mineral_species


def generate_solution():
    # List each species alongside their original values.
    species = original_values()[0]

    # Calculate the randomization for each specie.
    hco3 = np.random.uniform(1e-3, 2e-1)
    na = np.random.uniform(1e-2, 3)
    cl = na * species["Cl"] / species["Na"]
    ca = np.random.uniform(1e-4, 5e-3)
    mg = np.random.uniform(1e-5, 1e-3)
    s = np.random.uniform(1e-4, 5e-2)

    # Randomized sample.
    sample = {
        "{HCO3}": f"{hco3:.4e}",
        "{NA}": f"{na:.4e}",
        "{CL}": f"{cl:.4e}",
        "{CA}": f"{ca:.4e}",
        "{MG}": f"{mg:.4e}",
        "{S(6)}": f"{s:.4e}",
    }

    return sample


def generate_gas(volume, porosity):
    GAS_PRESSURE = 60  # bar

    volume_liquid = volume * porosity * 1000  # L
    fraction_water = np.random.uniform(0.15, 0.85)
    volume_water = volume_liquid * fraction_water  # L
    volume_gas = volume_liquid - volume_water  # L

    # List each species alongside their original values.
    species = original_values()[1]

    # Calculate the randomization for each specie.
    pressure_h2 = species["H2"]
    pressure_co2 = np.random.uniform(1e-6, 10)
    pressure_h2s = np.random.uniform(1e-6, 1e-2)
    pressure_n2 = np.random.uniform(1e-6, 10)
    pressure_ch4 = (
        GAS_PRESSURE - pressure_co2 - pressure_h2s - pressure_n2 - pressure_h2
    )

    # Randomized sample.
    sample = {
        "{CH4}": f"{pressure_ch4:.4e}",
        "{CO2}": f"{pressure_co2:.4e}",
        "{N2}": f"{pressure_n2:.4e}",
        "{H2S}": f"{pressure_h2s:.4e}",
    }

    return sample, volume_gas, volume_water


def generate_mineral(volume, porosity):
    # Fraction, Diameter (m), Density (kg/m^3), MW (kg/mol)
    species = original_values()[2]

    mass = np.array(list(species.values()))[:, 2] * volume * (1 - porosity)
    surface_area = (
        (6 / np.array(list(species.values()))[:, 1])
        * (mass / np.array(list(species.values()))[:, 2])
        * 1e-3
    )
    mol = mass / np.array(list(species.values()))[:, 3] * 1e-3

    sample_surface_area = {
        "{CALCITE_AREA}": f"{surface_area[0]:.4e}",
    }

    sample_mol = {
        "{CALCITE_MOLES}": f"{mol[0]:.4e}",
    }

    return sample_surface_area, sample_mol


def work(iter_num, model):
    print(f"Simulation {iter_num} is running...")

    porosity = np.random.uniform(0.1, 0.3)
    VOLUME = 1  # m3
    temp = np.random.uniform(25, 60)  # degC

    sample_solution = generate_solution()
    sample_gas, volume_gas, volume_water = generate_gas(VOLUME, porosity)
    (
        sample_mineral_surface_area,
        sample_mineral_mol,
    ) = generate_mineral(VOLUME, porosity)

    sample_props = {
        "{TEMPERATURE}": f"{temp:.4e}",
        "{WATER_VOLUME}": f"{volume_water:.4e}",
        "{GAS_VOLUME}": f"{volume_gas:.4e}",
    }

    sample_total = (
        sample_props
        | sample_solution
        | sample_gas
        | sample_mineral_mol
        | sample_mineral_surface_area
    )
    sample_list = [[k, v] for k, v in sample_total.items()]

    input_file = f"./input/{iter_num}_Input.txt"

    # Save the generated simulation results.
    np.savetxt(input_file, sample_list, fmt="%s")
    try:
        output, logs = run_model(model, sample_list)
        np.savetxt(f"./output/{iter_num}_Output.txt", output, fmt="%s")
        print(f"Simulation {iter_num} completed.")

        if DEBUG:
            generate_plots(iter_num, f"./output/{iter_num}_Output.txt")

    except subprocess.TimeoutExpired:
        logs = f"Simulation {iter_num} timed out after {MAX_DURATION} seconds."
        print(f"Simulation {iter_num} timed out.")

    except (RuntimeError, ValueError) as e:
        logs = str(e)
        print(f"Simulation {iter_num} failed.")

    if DEBUG:
        with open(f"./logs/{iter_num}_log.txt", "w") as f:
            f.write(logs)


def run_model(model, inputs):
    if sys.platform == "linux":
        DATABASE = r"/usr/share/phreeqc/database/phreeqc.dat"
        PHREEQC = r"phreeqc"
        with open(model, "r", encoding="latin-1") as file:
            raw = file.read()
    else:
        DATABASE = r"C:\phreeqc\database\phreeqc.dat"
        PHREEQC = r"C:\phreeqc\phreeqc.exe"
        with open(model, "r") as file:
            raw = file.read()

    # Inputs is a list of lists passed from work()
    for key, value in inputs:
        raw = raw.replace(key, str(value))

    with tempfile.TemporaryDirectory() as temp_dir:
        input = os.path.join(temp_dir, "model.pqi")
        output = os.path.join(temp_dir, "model.pqo")
        punch = os.path.join(temp_dir, "selected_output_1.sel")

        with open(input, "w") as f:
            f.write(raw)

        cmd = [PHREEQC, input, output, DATABASE]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=temp_dir,
            check=False,
            timeout=MAX_DURATION,
        )

        # Save execution logs
        with open(output, "r") as f:
            logs = f.read()

        if result.stderr:
            logs += f"\n\n[STDERR]:\n{result.stderr.strip()}"

        if result.returncode != 0:
            raise RuntimeError(logs)

        # Parse and return output
        output = []
        if os.path.exists(punch):
            with open(punch, "r") as f:
                output = [line.strip().split() for line in f]

    return output, logs


def generate_plots(iter_num, file_path):
    try:
        with open(file_path, "r") as f:
            headers = f.readline().strip().split()

        col_idx = {h: i for i, h in enumerate(headers)}

        data = np.genfromtxt(file_path, skip_header=1)
        time_d = data[1:, col_idx["time_d"]]

        fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(14, 15))
        fig.suptitle(
            f"Simulation Results - Iteration {iter_num}", fontsize=14, fontweight="bold"
        )

        # USER_GRAPH 1: pH
        axes[0, 0].plot(time_d, data[1:, col_idx["pH"]], color="blue")
        axes[0, 0].set_title("pH")
        axes[0, 0].set_xlabel("Time (days)")
        axes[0, 0].set_ylabel("pH")
        axes[0, 0].grid(True, linestyle="--", alpha=0.7)

        # USER_GRAPH 2: Trona-K model
        graph2_mapping = {
            "H2_g_mol": "H2(g)",
            "CO2_g_mol": "CO2(g)",
            "SO4_mol": "SO4",
            "Formate": "Formate",
            "Acetate": "Acetate",
            "Ca_mol": "Ca",
        }
        for col, label in graph2_mapping.items():
            if col in col_idx:
                axes[0, 1].plot(time_d, data[1:, col_idx[col]] * 1000, label=label)
        axes[0, 1].set_title("Trona-K model")
        axes[0, 1].set_xlabel("Time (days)")
        axes[0, 1].set_ylabel("Amount (mmol)")
        axes[0, 1].legend(loc="best", fontsize="small", ncol=2)
        axes[0, 1].grid(True, linestyle="--", alpha=0.7)

        # USER_GRAPH 3: Biomass growth
        graph3_species = ["X_SRB_mol", "X_IRB_mol"]
        for spec in graph3_species:
            if spec in col_idx:
                axes[1, 0].plot(time_d, data[1:, col_idx[spec]], label=spec)
        axes[1, 0].set_title("Biomass (moles of kinetic reactant)")
        axes[1, 0].set_xlabel("Time (days)")
        axes[1, 0].set_ylabel("Moles")
        axes[1, 0].legend(loc="best", fontsize="small")
        axes[1, 0].grid(True, linestyle="--", alpha=0.7)

        # USER_GRAPH 4: Ca Fe HS-
        graph4_mapping = {"Ca_mol": "Ca_model", "Fe_mol": "Fe_model", "HS-_mol": "HS-"}
        for col, label in graph4_mapping.items():
            if col in col_idx:
                axes[1, 1].plot(time_d, data[1:, col_idx[col]] * 1000, label=label)
        axes[1, 1].set_title("Ca Fe HS-")
        axes[1, 1].set_xlabel("Time (days)")
        axes[1, 1].set_ylabel("Amount (mmol)")
        axes[1, 1].legend(loc="best", fontsize="small")
        axes[1, 1].grid(True, linestyle="--", alpha=0.7)

        # USER_GRAPH 5: Iron reduction
        graph5_species = ["IR_H2", "IR_AC", "IR_FOR"]
        for spec in graph5_species:
            if spec in col_idx:
                axes[2, 0].plot(time_d, data[1:, col_idx[spec]], label=spec)
        axes[2, 0].set_title("Iron reduction")
        axes[2, 0].set_xlabel("Time (days)")
        axes[2, 0].set_ylabel("Kinetics")
        axes[2, 0].legend(loc="best", fontsize="small")
        axes[2, 0].grid(True, linestyle="--", alpha=0.7)

        # USER_GRAPH 6: Solution (Totals & HCO3-)
        graph6_mapping = {
            "HCO3_mol": "HCO3-",
            "Na_tot": "Na_tot",
            "Mg_tot": "Mg_tot",
            "Cl_tot": "Cl_tot",
            "Ca_tot": "Ca_tot",
            "S6_tot": "S6_tot",
        }
        for col, label in graph6_mapping.items():
            if col in col_idx:
                axes[2, 1].plot(time_d, data[1:, col_idx[col]] * 1000, label=label)
        axes[2, 1].set_title("Solution")
        axes[2, 1].set_xlabel("Time (days)")
        axes[2, 1].set_ylabel("Amount (mmol)")
        axes[2, 1].legend(loc="best", fontsize="small", ncol=2)
        axes[2, 1].grid(True, linestyle="--", alpha=0.7)

        plt.tight_layout(rect=(0, 0, 1, 0.96))

        plt.savefig(f"./plot/{iter_num}_Plot.jpg", dpi=300)

        plt.clf()
        plt.close(fig)

    except (OSError, ValueError, KeyError, IndexError) as e:
        print(f"Error plotting iteration {iter_num}: {e}")


def parse_inputs():
    input_dir = Path("./input")

    params = {}
    str_tracker = None

    # Gather all the values in inputs into a single dictionary.
    for file in input_dir.iterdir():
        contents = Path.read_text(file).split()

        for content in contents:
            try:
                content = float(content)

                if str_tracker is not None:
                    params[str_tracker].append(content)
            except ValueError:
                str_tracker = content

            if str_tracker not in params:
                params[str_tracker] = []

    # Parse the params dictionary.
    for key, value in params.items():
        params[key] = (float(np.min(value)), float(np.max(value)))

    return params


def generate_report(WORKERS, MODEL_SIZE):
    parsed_inputs = parse_inputs()

    output_count = sum(1 for item in Path("./output").iterdir() if item.is_file())

    with open("report.txt", "w") as f:
        separator_count = 80
        f.write(
            separator_count * "="
            + "\n"
            + "General Information"
            + "\n"
            + separator_count * "="
            + 2 * "\n"
            + f"Simulated Model => {MODEL}\n"
            + f"Thread Count => {WORKERS}\n"
            + f"Simulation Timeout => {MAX_DURATION} seconds\n"
            + f"Debug Flag => {DEBUG}\n"
            + f"Number of Simulations => {MODEL_SIZE}\n"
            + f"Number of Successful Simulations => {output_count}\n"
            + f"Number of Failed Simulations => {MODEL_SIZE - output_count}\n"
            + "\n"
            + separator_count * "="
            + "\n"
            + "Input Parameter Ranges"
            + "\n"
            + separator_count * "="
            + 2 * "\n"
        )
        for key, value in parsed_inputs.items():
            f.writelines([str(key), " => ", str(value), "\n"])

    print("Report have been generated.")


def organize_folders():
    input_dir = r"./input"
    output_dir = r"./output"
    failed_dir = r"./failed"

    os.makedirs(failed_dir, exist_ok=True)

    inputs = {
        f.split("_")[0] for f in os.listdir(input_dir) if f.endswith("_Input.txt")
    }
    outputs = {
        f.split("_")[0] for f in os.listdir(output_dir) if f.endswith("_Output.txt")
    }

    missing = inputs - outputs

    for m in sorted(missing, key=int):
        file_name = f"{m}_Input.txt"
        src_path = os.path.join(input_dir, file_name)
        dst_path = os.path.join(failed_dir, file_name)

        shutil.move(src_path, dst_path)
        print(f"Moved to failed: {file_name}")


if __name__ == "__main__":
    main()
