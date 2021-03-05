import os
import numpy as np
from pyscf import gto, scf, mcscf
from pyscf.shciscf import shci
from pyscf.tools import molden
from pyscf.lib import chkfile

#
# User settings
#
eps1 = 1e-3
ncas, nelecas = (100, 100)
outputbase = f"1_C_{nelecas}e_{ncas}o_eps1={eps1}"
uhf_outputbase = "uhf_1_C"
scratchdir = os.environ.get("TMP")
nprocs = os.environ.get("OMP_NUM_THREADS")  # This is used for MPI in Dice

# Set up data directories
data_dirs = ["_logs", "_chk", "_molden", "_data"]
for dd in data_dirs:
    if not os.path.exists(dd):
        os.mkdir(dd)


# Build Molecule
mol = gto.Mole()
mol.atom = """Fe 0.025745 -0.124890 0.221534
N 1.973404 0.308372 -0.012535
N -0.016886 1.852833 0.129448
N -1.992657 0.322677 -0.234778
N -0.374202 -1.341364 1.594130
N -0.617028 -2.095226 2.391073
C 2.298988 1.606163 0.028942
C 1.187405 2.499698 0.029803
C -1.182193 2.498217 -0.133265
C -2.305587 1.603075 -0.251088
C 2.984541 -0.678749 -0.060024
C 3.339409 -1.222013 -1.303465
C 4.291921 -2.233778 -1.339411
H 4.570862 -2.656491 -2.297734
C 4.879050 -2.702602 -0.176176
H 5.616639 -3.493367 -0.220286
C 4.513549 -2.160762 1.044823
C 3.566096 -1.146606 1.128112
C -2.953994 -0.703887 -0.272069
C -2.841034 -1.651663 -1.304572
C -3.737413 -2.708940 -1.346104
C -4.708667 -2.858000 -0.368352
H -5.397627 -3.691452 -0.408097
C -4.770128 -1.954222 0.676135
H -5.497187 -2.092300 1.468056
C -3.900453 -0.869651 0.752425
H -3.660569 -3.429868 -2.151492
H 4.964596 -2.529790 1.958675
C 2.700009 -0.711454 -2.553965
H 3.022428 -1.280615 -3.423377
H 1.609935 -0.760663 -2.495519
H 2.938298 0.338653 -2.734277
C 3.163821 -0.569843 2.446302
H 3.452357 0.479170 2.540384
H 2.081319 -0.595778 2.579620
H 3.623608 -1.112672 3.269261
C -1.744161 -1.527410 -2.310395
H -0.760108 -1.646878 -1.836896
H -1.821867 -2.292411 -3.080251
H -1.727116 -0.550371 -2.794101
C -3.952303 0.045866 1.933921
H -4.288933 -0.491318 2.818809
H -2.977887 0.481084 2.154750
H -4.647895 0.874580 1.785483
C 3.715242 2.063833 0.003418
H 3.792051 3.144841 0.077712
H 4.292937 1.631881 0.822954
H 4.218109 1.753994 -0.916005
C -3.682021 2.137366 -0.451622
H -4.042367 2.665879 0.433016
H -3.694396 2.857238 -1.270815
H -4.389906 1.346264 -0.685104
C 1.201904 3.885041 -0.133510
C -1.200471 3.873817 -0.316113
H 2.142147 4.418723 -0.168916
H -2.132840 4.391273 -0.497194
C 0.006644 4.571752 -0.281280
H 0.014687 5.645010 -0.414719
"""
mol.basis = "ccpvdz"
mol.verbose = 6
mol.output = f"_logs/{outputbase}.out"
mol.max_memory = 100000  # in MB
mol.build()


mf = scf.UHF(mol).newton()
mf.__dict__.update(chkfile.load(f"_chk/{uhf_outputbase}.chk", "scf"))

#
# HCISCF Calculation
#
mc = shci.SHCISCF(mf, ncas, nelecas).density_fit()
mc.chkfile = f"_chk/{outputbase}.chk"
mc.natorb = True
mc.fcisolver.sweep_iter = [0, 3, 6, 9]
mc.fcisolver.sweep_epsilon = [eps1] * 4
mc.fcisolver.mpiprefix = f"mpirun -np {nprocs}"
mc.fcisolver.scratchDirectory = f"{scratchdir}/fepdi/singlet/{outputbase}"
if not os.path.exists(mc.fcisolver.scratchDirectory):
    # makedirs creates subfolders too
    os.makedirs(mc.fcisolver.scratchDirectory)

# Use UNO orbitals as initial guess
uhf_NOs = np.load(f"_chk/{uhf_outputbase}_NOs.npy")
mc.conv_tol = 1e-4
mc.mc1step(uhf_NOs)

# Save Data
np.savetxt(f"_data/{outputbase}_noons.txt", mc.mo_occ)
molden.from_mo(mol, f"_molden/{outputbase}.molden", mc.mo_coeff)

# Print population info
mc.analyze()
