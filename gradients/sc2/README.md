# Sc Dimer Gradients as a Function of Method Flavor

Calculations to generate the $Sc_2$ Single Pt. figure where we show the error in gradients as a function of which method we use.

## Reproducing the Data

```bash
cd single_point
sh run_all.sh
cd ../eps1_convergence
sh run_all.sh
```

## Results

- FCI space size is 154224 determinants (setting $\epsilon_1 = 0$)

![](single_point_error.png)

- The results vary somewhat from run to run for the active-active rotations, in particular the hciscf_aa error varies in the first sig fig.
  - Below are output of `plot_grad.py` from two different calculations of vhciscf_aa and hciscf_aa.
