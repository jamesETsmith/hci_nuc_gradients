# Gradients and Geometry Optimization with HCISCF Wave Functions

This repository contain all files required to reproduce the data reported in the paper title "TITLE". [ADD CITATION]

If this work helped your own research efforts, please cite our work using:

```
Citation coming soon!
```

## General Organization

```
.
├── geometry_opt
│   ├── dft
│   └── hciscf
├── gradients
│   ├── finite_difference
│   ├── n2
│   └── sc2
└── README.md
```


## How to use this project
Each of the directories withing `./gradients` and `./geometry_opt` contain `README.md` files with detailed instructions to reproduce data from this project.
In general, most directories contain bash scripts called `run_all.sh` which automate the lowest level tasks.


## Software Version Details and Requirements

### DFT 
- Gaussian16
### HCISCF
- PySCF

### Analysis
- NumPy
- Matplotlib
- Seaborn
- Pandas
- rmsd
- cclib


