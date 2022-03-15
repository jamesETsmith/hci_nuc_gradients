# How To Use This Section

1. Run UHF survey
2. Select lowest energy UHF solutions
3. Run (100e,100o) UNO-HCISCF
4. Check UNO calcs (see `analysis/dz_check_uno_calcs`)
5. Launch geometry optimizations
6. Check the observables (Fe charge & NOONs) at the initial geometry
7. Visualize orbitals at the initial geometry
8. If 6 and 7 look good, run extrapolation to improve energy estimate


The majority of data from the paper comes from `dz` where we use the `cc-pVDZ` basis set.
During our investigation, we attempted to use `cc-pVTZ`, but often ran into hardware limitations.
Those files are included in `tz`.