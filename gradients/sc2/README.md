# Sc$_2$ Gradients Data

Calculations to generate the $Sc_2$ Single Pt. figure where we show the error in gradients as a function of which method we use.

## Reproducing the Data

```bash
cd single_point
sh run_all.sh
cd ../eps1_convergence
sh run_all.sh
```
> NOTE: If you want to reproduce the $\epsilon_1$ convergence data for N$_2$, you should run [`single_point`](single_point/) first (using [`run_all.sh`](single_point/run_all.sh)) because the [plotting script](eps1_convergence/plot_conv.py) relies on some of the data produced during those calculations.