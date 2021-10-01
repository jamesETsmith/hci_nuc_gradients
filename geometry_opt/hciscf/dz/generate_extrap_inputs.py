"""
Generates input.dat files for Dice calculations.
"""
import os
import numpy as np

singlet_occ = """"0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 """
triplet_occ = """"0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 """

input_header = """
#system
nocc {}
{} 
end
orbitals ../{}
nroots 1

#variational
schedule 
"""

input_footer = """

davidsonTol 5e-05
dE 1e-08
maxiter 21

#pt
nPTiter 1000
epsilon2 5.0e-08
epsilon2Large 1e-5
targetError 1e-5
sampleN 500

#misc
noio
DoRDM
DoSpinRDM
"""


def write_schedule(eps1: float) -> str:
    schedule = ""
    for i in range(10):
        schedule += f"{i}   {(9-i)*1e-5 + eps1:.2e}\n"
    return schedule + "end\n\n"


def write_input(
    eps1s: np.ndarray, ncas: int, spin: int, initial: bool, dice_path: str, n_proc: int
):

    # Choose spin
    if spin == 0:
        occ = singlet_occ
    elif spin == 2:
        occ = triplet_occ

    # Set FCIDUMP file name
    if initial:
        fcidump = "FCIDUMP_INITIAL"
        output_dir = "initial_extrapolation"
    else:
        fcidump = "FCIDUMP_FINAL"
        output_dir = "final_extrapolation"

    os.makedirs(output_dir, exist_ok=True)
    write_run_all(output_dir, dice_path, n_proc, eps1s.size)

    # Write input files for each epsilon
    for i, eps1 in enumerate(eps1s):
        with open(f"{output_dir}/input_{i}.dat", "w") as f:
            f.write(input_header.format(ncas, occ, fcidump))
            f.write(write_schedule(eps1))
            f.write(input_footer)


def write_run_all(output_dir: str, dice_path: str, n_proc: int, n_eps1: int):
    header = f"""#!/usr/bin/bash

DICE_EXE={dice_path}

# Run all Dice jobs
"""
    with open(f"{output_dir}/run_all.sh", "w") as f:
        f.write(header)
        # Run the most expensive ones first
        for i in reversed(range(n_eps1)):
            f.write(
                f"""mpirun -np {n_proc} $DICE_EXE input_{i}.dat > output_{i}.dat\n"""
            )
        f.write("\n")


if __name__ == "__main__":
    import argparse

    # fmt: off
    parser = argparse.ArgumentParser(description="""This script generates Dice input files for a set of epsilon_1 values so we can extrapolate to the FCI-limit.""")
    parser.add_argument("--initial", help="Initial geometry or final", type=int, choices=[1, 0], required=True)
    parser.add_argument("--ncas", help="The number of active space orbitals, which will be the same number as the electons in the CAS calc.", type=int, required=True)
    parser.add_argument("--spin", help="N_alpha - N_beta electrons.", type=int, choices=[0, 2], required=True)
    parser.add_argument("--dice_path", help="Path to the Dice executable", type=str, required=True)
    parser.add_argument("--n_proc", help="Number of MPI processes", type=int, required=True)
    args = parser.parse_args()
    # fmt: on

    # Shorthand for some vars
    ncas = args.ncas
    spin = args.spin
    initial = bool(args.initial)
    dice_path = args.dice_path
    n_proc = args.n_proc
    print(
        f"Generating Dice inputs for ncas = {ncas} spin={spin} initial geometry = {initial}"
    )
    print(f"Using Dice at {dice_path} with {n_proc} mpi processes")

    eps1s = np.linspace(5e-5, 1e-5, num=10)
    print("Epsilon_1 values", eps1s)

    write_input(eps1s, ncas, spin, initial, dice_path, n_proc)
