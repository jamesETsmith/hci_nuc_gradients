# UHF Survey

This directory contains files necessary to reproduce the survey of SCF options used to find the lowest energy UHF solution for use with multireference calculations.
This survey uses the `cc-pVTZ` basis set (as opposed to `cc-pVDZ`).

## Reproducing the Data

Users should modify `submit.sh` to fit their cluster scheduler.

```
./setup_dirs.sh
./launch_all.sh
```