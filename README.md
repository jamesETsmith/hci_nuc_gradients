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
- Matplotlib
- Seaborn
- Pandas
- rmsd
- cclib


## TODOs

### Gradients
- [X] FD Tests
- [ ] FD as a function of epsilon_1
- [X] N2 Single point
- [X] N2 eps1 convergence
- [X] Sc2 Single point
- [X] Sc2 eps1 convergence

### DFT
- [X] M06-L 
- [X] M06-2X
- [X] MN15  


### HCISCF
- [X] Fe(PDI) Big CAS (1,A)
- [ ] Fe(PDI) Big CAS (3,B)
- [X] Fe(PDI) Big CAS (3,C)

- Fe(PDI) HCISCF (1,A)
  - [X] (10e,10o) (Running)
  - [X] (20e,20o) (Running)
  - [X] (30e,30o) (Running)
  - [X] (40e,40o) (Running)
- Fe(PDI) HCISCF (3,A)
  - [X] (10e,10o) (Running)
  - [X] (20e,20o) (Running)
  - [X] (30e,30o) (Running)
  - [ ] (40e,40o) (Running)
- Fe(PDI) HCISCF (1,C)
  - [X] (10e,10o)
  - [X] (20e,20o)
  - [ ] (30e,30o)
  - [X] (40e,40o)