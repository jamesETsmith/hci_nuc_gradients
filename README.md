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

- PySCF
- Matplotlib
- Seaborn
- Pandas


## TODOs

### Gradients
- [X] FD Tests
- [ ] N2 Single point
- [ ] N2 eps1 convergence
- [ ] Sc2 Single point
- [ ] Sc2 eps1 convergence

### DFT
- [ ] M06-L
- [ ] M06-2X
- [ ] MN15


### HCISCF
- [ ] Fe(PDI) Big CAS (SSG)
- [ ] Fe(PDI) Big CAS (TSG)

- Fe(PDI) HCISCF (1,SSG)
  - [ ] (10e,10o)
  - [ ] (20e,20o)
  - [ ] (30e,30o)
  - [ ] (40e,40o)
- Fe(PDI) HCISCF (3,SSG)
  - [ ] (10e,10o)
  - [ ] (20e,20o)
  - [ ] (30e,30o)
  - [ ] (40e,40o)
- Fe(PDI) HCISCF (1,TSG)
  - [ ] (10e,10o)
  - [ ] (20e,20o)
  - [ ] (30e,30o)
  - [ ] (40e,40o)