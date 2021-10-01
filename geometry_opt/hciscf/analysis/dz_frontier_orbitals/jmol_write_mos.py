"""
Wrote PNG images of MOs rendered with JMOL.
"""
import os
import subprocess


def run_jmol_from_molden(
    molden_file: str, java_path: str, jmol_path: str, output_dir: str
):

    rotate_molecule = ""
    if "3_C" in molden_file or "TSG" in molden_file:
        rotate_molecule = "moveto 0 AXIS c"

    jmol_script_template = f"""initialize;
set background [xffffff];
load {molden_file}
{rotate_molecule}

# Set MO display settings
mo fill
mo nomesh
mo cutoff 0.02
mo TITLEFORMAT "%F  MO=%I"

# Load and write MOs
mo  118; write IMAGE 600 600 PNG 0 "{output_dir}/118.png"
mo  119; write IMAGE 600 600 PNG 0 "{output_dir}/119.png"
mo  120; write IMAGE 600 600 PNG 0 "{output_dir}/120.png"
mo  121; write IMAGE 600 600 PNG 0 "{output_dir}/121.png"

exitJmol

    """

    script_name = "jmol_script.spt"
    with open(script_name, "w") as f:
        f.write(jmol_script_template)

    subprocess.run([java_path, "-jar", jmol_path, script_name])
    os.remove(script_name)


def run_jmol_from_cube(functional: str, species: str, java_path: str, jmol_path: str):

    species_base = functional + "_" + species
    rotate_molecule = ""
    if species == "3_C":
        rotate_molecule = "moveto 0 AXIS c"

    jmol_script_template = f"""initialize;
set background [xffffff];
load "_logs/{species_base}.log"
rotate z 180
{rotate_molecule}
isosurface sign red blue "_cubes/{species_base}_118.cub"; write IMAGE 600 600 PNG 0 "_figures/{species}/{functional}/118.png"
isosurface sign red blue "_cubes/{species_base}_119.cub"; write IMAGE 600 600 PNG 0 "_figures/{species}/{functional}/119.png"
isosurface sign red blue "_cubes/{species_base}_120.cub"; write IMAGE 600 600 PNG 0 "_figures/{species}/{functional}/120.png"
isosurface sign red blue "_cubes/{species_base}_121.cub"; write IMAGE 600 600 PNG 0 "_figures/{species}/{functional}/121.png"

exitJmol
"""

    # Make sure the _figures directory is there
    os.makedirs(os.path.join("_figures", species, functional), exist_ok=True)

    # Write script
    script_name = "jmol_script.spt"
    with open(script_name, "w") as f:
        f.write(jmol_script_template)

    # Run script
    subprocess.run([java_path, "-jar", jmol_path, script_name])
    os.remove(script_name)


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 5:
        raise AssertionError(
            "Wrong number of command line args. Try\n"
            + "python jmol_write_mos.py my_orbs.molden java_path jmol_jar_path output_dir"
        )

    run_jmol_from_molden(*sys.argv[1:])
