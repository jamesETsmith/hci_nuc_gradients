import os
import numpy as np


def parse_for_noons(filename: str):

    STOP_TOK = "Alpha Molecular Orbital Coefficients:"
    noons = {"Alpha": [], "Beta": []}

    # Read the file and grab the NOONs for Alpha and Beta
    with open(filename, "r") as f:
        done_reading = False
        while not done_reading:
            line = f.readline()

            if STOP_TOK in line:
                done_reading = True
                break

            if "occ. eigenvalues" in line or "virt. eigenvalues" in line:
                # print(line.split())
                lsplit = line.split()
                noons[lsplit[0]] += list(map(float, lsplit[4:]))

    # Check that the NOONs for Alpha and Beta are equal
    noons["Alpha"] = np.array(noons["Alpha"])
    noons["Beta"] = np.array(noons["Beta"])
    diff = np.linalg.norm(noons["Alpha"] - noons["Beta"])
    # print(noons)
    if diff > 1e-14:
        raise AssertionError(
            "Alpha and Beta NOONs don't match " + f"l-2 norm of the difference = {diff}"
        )

    # Save the NOONs
    os.makedirs("_data", exist_ok=True)
    output_file = "_data/" + os.getcwd().split("/")[-1] + "_noons.txt"
    np.savetxt(output_file, noons["Alpha"])


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        raise ValueError(
            "Wrong number of CLIs. Try python parse_for_noons.py my_log_file.log"
        )

    parse_for_noons(*sys.argv[1:])
