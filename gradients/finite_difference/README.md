# N2 Gradient and Finite Differece (FD) Tests

## Directory Manifest

The subdirectories (`casscf`, `vhciscf`, and `hciscf`) contain:
- `n2_casscf_grad.py` which runs a single FD calculation for a specified FD step size
- `run_all.sh` runs `n2_casscf_grad.py` for a set of step size and performs a simple grep command to test for convergence.
- `plot_fd_grad_conv.py` Reads the `*.json` outputs and generates a convergence plot.
- `fd_comparison.png` Shows the difference between analytical gradients and FD ones as a function of step size


## Reproducing Results
- In each of the subdirectories (`casscf`, `vhciscf`, and `hciscf`), run `run_all.sh`. 
- (Optional) run `plot_fd_grad_conv.py` and visualize the convergence of error w.r.t. the finite difference step size.
- In the current directory, run `compare_fd_error.py` to collect a subset of the data and visually compare them. The outputs are saved in `_figures`.


