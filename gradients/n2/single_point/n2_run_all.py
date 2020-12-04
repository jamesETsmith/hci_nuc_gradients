"""
This script calculates the N2 gradients using CASSCF.
See `README.md` for more details.
"""
import numpy as np
from pyscf.shciscf import shci

import util

# Control which calculations to run
CASSCF = 1
vHCISCF = 1
vHCISCF_AA = 1
HCISCF = 1
HCISCF_AA = 1

# Set HCI epsilon1 parameter
epsilon1 = 5e-3

#
# CASSCF
#
if CASSCF:
    mc = util.make_n2("_logs/casscf.out")
    mc.mc2step()

    grad = mc.Gradients().kernel()
    np.save("_data/grad_casscf", grad)
    del mc

#
# vHCISCF
#
if vHCISCF:
    mc2 = util.make_n2("_logs/vhciscf.out")
    mc2.fcisolver = shci.SHCI(mc2.mol)
    mc2.fcisolver.sweep_epsilon = [epsilon1]
    mc2.fcisolver.sweep_iter = [0]
    mc2.mc2step()

    grad2 = mc2.Gradients().kernel()
    np.save("_data/grad_vhciscf", grad2)
    del mc2

#
# vHCISCF + AA
#
if vHCISCF_AA:
    mc3 = util.make_n2("_logs/vhciscf_aa.out")
    mc3.fcisolver = shci.SHCI(mc3.mol)
    mc3.fcisolver.sweep_epsilon = [epsilon1]
    mc3.fcisolver.sweep_iter = [0]
    mc3.mc2step()

    mc3.internal_rotation = True
    mc3.max_cycle_macro = 300
    mc3.conv_tol = 1e-6
    mc3.conv_tol_grad = np.sqrt(mc3.conv_tol)
    mc3.mc2step()

    grad3 = mc3.Gradients().kernel()
    np.save("_data/grad_vhciscf_aa", grad3)

    del mc3


#
# HCISCF
#
if HCISCF:
    mc4 = util.make_n2("_logs/hciscf.out")
    mc4.fcisolver = shci.SHCI(mc4.mol)
    mc4.fcisolver.sweep_epsilon = [epsilon1]
    mc4.fcisolver.sweep_iter = [0]
    mc4.fcisolver.stochastic = False
    mc4.fcisolver.epsilon2 = 1e-10
    mc4.mc2step()

    grad4 = mc4.Gradients().kernel()
    np.save("_data/grad_hciscf", grad4)

    del mc4
#
# HCISCF + AA
#
if HCISCF_AA:
    mc5 = util.make_n2("_logs/hciscf_aa.out")
    mc5.fcisolver = shci.SHCI(mc5.mol)
    mc5.fcisolver.sweep_epsilon = [epsilon1]
    mc5.fcisolver.sweep_iter = [0]
    mc5.fcisolver.stochastic = False
    mc5.fcisolver.epsilon2 = 1e-10
    mc5.mc2step()

    mc5.internal_rotation = True
    mc5.max_cycle_macro = 300
    mc5.conv_tol = 3e-6
    mc5.conv_tol_grad = np.sqrt(mc5.conv_tol)
    mc5.mc2step()

    grad5 = mc5.Gradients().kernel()
    np.save("_data/grad_hciscf_aa", grad5)

    del mc5
