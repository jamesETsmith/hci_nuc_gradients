import os
import numpy as np
from pyscf import gto, scf, dft
from pyscf.tools import molden
from pyscf.lib import chkfile
from pyscf.scf.uhf import mulliken_spin_pop
from scipy.linalg import eigh
from functools import reduce

#
# User Settings
#
outputbase = "uhf_1_B"


# Set up data directories
data_dirs = ["_logs", "_chk", "_molden", "_data"]
for dd in data_dirs:
    if not os.path.exists(dd):
        os.mkdir(dd)


# Build Molecule
mol = gto.Mole()
mol.atom = """Fe   -0.000034    0.022319    0.000241
N          1.907328    0.362588    0.002684
N          0.000004    1.891890    0.000191
N         -1.907303    0.362612   -0.002062
N         -0.000005   -1.789275    0.000030
N         -0.000140   -2.911223   -0.000153
C          2.289904    1.641764   -0.000498
C          1.198098    2.556878   -0.000328
C         -1.198084    2.556898    0.000721
C         -2.289906    1.641795    0.000852
C          2.880727   -0.668393   -0.000580
C          3.306166   -1.198975   -1.223414
C          4.213571   -2.252120   -1.204901
H          4.548790   -2.671678   -2.146131
C          4.678239   -2.771871   -0.008242
H          5.382166   -3.593716   -0.011183
C          4.226134   -2.249930    1.192241
C          3.318851   -1.196764    1.218588
C         -2.880728   -0.668379    0.000595
C         -3.319004   -1.195926   -1.218908
C         -4.226270   -2.249101   -1.193143
C         -4.678221   -2.771865    0.007055
H         -5.382134   -3.593723    0.009500
C         -4.213410   -2.252933    1.203997
H         -4.548504   -2.673102    2.144996
C         -3.305998   -1.199784    1.223097
H         -4.571447   -2.666403   -2.131777
H          4.571209   -2.667817    2.130649
C          2.757424   -0.658426   -2.503353
H          3.140410   -1.205774   -3.361717
H          1.665120   -0.721445   -2.517924
H          2.997983    0.397141   -2.643256
C          2.783460   -0.654635    2.503485
H          3.023287    0.401697    2.638597
H          1.691727   -0.720188    2.530469
H          3.177480   -1.199379    3.358504
C         -2.783588   -0.653113   -2.503504
H         -1.691843   -0.718533   -2.530408
H         -3.177469   -1.197514   -3.358806
H         -3.023511    0.403247   -2.638209
C         -2.757225   -0.659969    2.503338
H         -3.140575   -1.207472    3.361439
H         -1.664949   -0.723452    2.518088
H         -2.997373    0.395652    2.643522
C          3.718403    2.051059   -0.004148
H          3.822623    3.132536   -0.009163
H          4.248752    1.666053    0.869692
H          4.246170    1.658122   -0.875993
C         -3.718396    2.051100    0.004866
H         -4.245578    1.659089    0.877499
H         -3.822597    3.132585    0.008836
H         -4.249361    1.665198   -0.868188
C          1.207972    3.946622   -0.000455
C         -1.207938    3.946635    0.000752
H          2.146814    4.483402   -0.000797
H         -2.146769    4.483434    0.001163
C          0.000026    4.637267    0.000127
H          0.000031    5.718597    0.000070
"""
mol.basis = "ccpvdz"
mol.verbose = 5
mol.output = f"_logs/{outputbase}.out"
mol.max_memory = 100000  # in MB
mol.build()


# Without restarting from UKS we get wrong spin density on
# Mean Field
mf0 = dft.UKS(mol, xc="PBE0")
mf0.diis = scf.ADIIS()
mf0.level_shift = 0.0001  # NEEDS TO BE SET TO SOMETHING SMALL
mf0 = scf.fast_newton(mf0)


mf = scf.UHF(mol).newton()
mf.diis = scf.ADIIS()
mf.chkfile = f"_chk/{outputbase}.chk"
mf.kernel(mf0.make_rdm1())

new_mo, _ = mf.stability()

# Check for stability

stable = False
if (np.array(new_mo) == np.array(mf.mo_coeff)).all():
    stable = True
else:
    dm = mf.make_rdm1(new_mo)
    mf.kernel(dm0=dm)

    new_mo, _ = mf.stability()
    if (np.array(new_mo) == np.array(mf.mo_coeff)).all():
        stable = True

mf.analyze()
scf.uhf.mulliken_spin_pop(mol, mf.make_rdm1(), s=mf.get_ovlp())

# Dump orbitals
molden.from_scf(mf, f"_molden/{outputbase}.molden")

# Dump objects for making NO
dm_a, dm_b = mf.make_rdm1()
S = mf.get_ovlp()

Dm = dm_a + dm_b
A = reduce(np.dot, (S, Dm, S))
w, v = eigh(A, b=S)
noons = np.flip(w)
NOs = np.flip(v, axis=1)

np.savetxt(f"_data/{outputbase}_noons.txt", noons)
np.save(f"_chk/{outputbase}_NOs.npy", NOs)
molden.from_mo(mol, f"_molden/{outputbase}_NOs.molden", NOs)
