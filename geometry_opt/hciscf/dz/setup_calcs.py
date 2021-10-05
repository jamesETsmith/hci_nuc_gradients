import os

for s in [0, 2]:
    m = 2 * (s // 2) + 1
    for g in ["A", "B", "C"]:
        for ncas in [10, 20, 30, 40, 50, 60]:
            os.makedirs(f"{m}_{g}/{ncas}e_{ncas}o", exist_ok=True)
            eps1 = 7.5e-5
            # if ncas == 40:
            #     eps1 = 2.0e-4
            if ncas >= 50:
                eps1 = 1.0e-4

            with open(f"{m}_{g}/{ncas}e_{ncas}o/run_all.sh", "w") as f:
                f.write(
                    f"""#!/usr/bin/bash

# Set up directories
mkdir -p _molden _logs _data _chk

# Geometry optimization
python ../../opt.py --eps1={eps1} --ncas={ncas} --spin={s} --geometry={g} 2>_logs/geometric_{eps1}.out
"""
                )
