# vHCISCF Geometry Optimization with cc-pVDZ
 
This directory contains files necessary to reproduce the generation of the UNO orbitals and the geometry optimizations where they're used.

## Coarse Manifest

- `opt.py`: The workhorse script for the geometry optimization. It takes several command line arguments and reads in UNO orbitals.
- `setup_calcs.py`: Generates bash scripts to run `opt.py` in the child directories of the current directory.
- `generate_extrap_inputs.py`: Generates the SHCI extrapolation input files at the initial and final geometries from the vHCISCF geometry optimization.

## Reproducing the Data

Users will need to run `python vhci.py` in all of the `M_X/gen_uno` directories in the current directory where M ∈ [A, B, C] and X ∈ [1, 3].
This will generate the UNO orbitals for the multireference calculations and requires that users have already run the necessary calculations in `../uhf_survey`.

Next users will need to generate the bash scripts necessary to run `opt.py`:

```
python setup_calcs.py
```

After generating (and running) those calculations, users will need to generate the Dice inputs for generating the SHCI extrapolation.
From every grandchild directory (e.g. `1_A/40e_40o`), users will need to run `generate_extrap_inputs.py`.
The command line execution should look something like:

```
python3 ../../generate_extrap_inputs.py --initial=1 --ncas=40 --spin=0 --dice_path=/mnt/home/jsmith/apps/pull_requests/Dice/Dice --n_proc=128
```