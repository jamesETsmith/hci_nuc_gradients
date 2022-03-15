# Gradients and Geometry Optimization with HCISCF Wave Functions

This repository contain all files required to reproduce the data reported in the paper ["Nuclear Gradients of Near-Exact Complete Active Space Self-Consistent Field Wave Functions"](https://arxiv.org/abs/2201.06514)

If this work helped your own research efforts, please cite our work using:

```
@article{smith2022nuclear,
      title={Nuclear Gradients of Near-Exact Complete Active Space Self-Consistent Field Wave Functions}, 
      author={James E. T. Smith and Joonho Lee and Sandeep Sharma},
      year={2022},
      eprint={2201.06514},
      archivePrefix={arXiv},
      primaryClass={physics.chem-ph}
}
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


