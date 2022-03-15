import os
import numpy as np
from pyscf import gto, scf, mcscf
from pyscf.shciscf import shci
from pyscf.tools import molden
from pyscf.lib import chkfile

#
# User settings
#
eps1 = 1e-4
ncas, nelecas = (100, 100)
outputbase = f"1_A_DZ_{nelecas}e_{ncas}o_eps1={eps1}_vhci"
uhf_outputbase = "../../../uhf_survey/FN_ADIIS/1_A/_chk/uhf_1_A_newton"
scratchdir = os.environ.get("TMP")
nprocs = os.environ.get("OMP_NUM_THREADS")  # This is used for MPI in Dice

# Set up data directories
data_dirs = ["_logs", "_chk", "_molden", "_data"]
for dd in data_dirs:
    os.makedirs(dd, exist_ok=True)


# Build Molecule
mol = gto.Mole()
mol.atom = """Fe -0.000000 0.017939 -0.000012
N 1.905716 0.352827 -0.000010
N 0.000000 1.902969 -0.000021 
N -1.905716 0.352828 0.000005
N 0.000000 -1.787363 0.000003
N -0.000003 -2.911341 0.000001
C 2.288283 1.636796 -0.000124 
C 1.197138 2.557818 -0.000086 
C -1.197137 2.557819 0.000040
C -2.288283 1.636796 0.000090
C 2.886288 -0.672389 0.000238 
C 3.325785 -1.197247 -1.220175 
C 4.244247 -2.240752 -1.197693 
H 4.589636 -2.655045 -2.137671 
C 4.707741 -2.757062 0.000710 
H 5.420562 -3.571363 0.000895 
C 4.244157 -2.240285 1.198875 
C 3.325706 -1.196756 1.220889 
C -2.886288 -0.672389 -0.000206 
C -3.325779 -1.196717 -1.220845 
C -4.244231 -2.240244 -1.198811 
C -4.707741 -2.757061 -0.000636 
H -5.420562 -3.571362 -0.000805 
C -4.244173 -2.240792 1.197756 
H -4.589499 -2.655123 2.137741 
C -3.325713 -1.197285 1.220218 
H -4.589588 -2.654159 -2.138968 
H 4.589450 -2.654238 2.139038 
C 2.782012 -0.662758 -2.504853 
H 3.203953 -1.184806 -3.360850 
H 1.694432 -0.770872 -2.544595 
H 2.983126 0.402990 -2.628223 
C 2.781670 -0.661965 2.505335 
H 2.982487 0.403865 2.628424 
H 1.694113 -0.770353 2.544992 
H 3.203638 -1.183645 3.361544 
C -2.781846 -0.661854 -2.505305 
H -1.694291 -0.770235 -2.545053 
H -3.203880 -1.183486 -3.361510 
H -2.982677 0.403984 -2.628313
C -2.781835 -0.662869 2.504883 
H -3.203710 -1.184965 3.360884 
H -1.694253 -0.770990 2.544534 
H -2.982934 0.402872 2.628335
C 3.716831 2.046181 -0.000083 
H 3.820546 3.127864 -0.003002 
H 4.245787 1.660251 0.874362 
H 4.247151 1.655263 -0.871439 
C -3.716830 2.046182 -0.000009
H -4.247179 1.655309 0.871349
H -3.820545 3.127865 0.002852
H -4.245758 1.660207 -0.874452
C 1.209585 3.949572 -0.000109
C -1.209584 3.949572 0.000035
H 2.147310 4.487836 -0.000142 
H -2.147309 4.487836 0.000057
C 0.000000 4.636142 -0.000043 
H 0.000000 5.717916 -0.000055
"""
mol.basis = "ccpvdz"
mol.verbose = 6
mol.output = f"_logs/{outputbase}.out"
mol.max_memory = 400000  # in MB
mol.build()


mf = scf.UHF(mol).newton()
mf.__dict__.update(chkfile.load(f"{uhf_outputbase}.chk", "scf"))

#
# vHCI Calculation
#
mc = mcscf.CASCI(mf, ncas, nelecas).density_fit()
mc.chkfile = f"_chk/{outputbase}.chk"
mc.natorb = True

mc.fcisolver = shci.SHCI(mol)
mc.fcisolver.sweep_iter = [0, 3, 6, 9]
mc.fcisolver.sweep_epsilon = [eps1] * 4
mc.fcisolver.mpiprefix = f"mpirun -np {nprocs}"
mc.fcisolver.scratchDirectory = f"{scratchdir}/fepdi/{outputbase}"
os.makedirs(mc.fcisolver.scratchDirectory, exist_ok=True)

# Use UNO orbitals as initial guess
uhf_NOs = np.load(f"{uhf_outputbase}_NOs.npy")
mc.kernel(uhf_NOs)

# Save Data
np.savetxt(f"_data/{outputbase}_noons.txt", mc.mo_occ)
np.save(f"_chk/{outputbase}_NOs.npy", mc.mo_coeff)
molden.from_mo(mol, f"_molden/{outputbase}.molden", mc.mo_coeff)

# Print poputlation info
mc.analyze()
