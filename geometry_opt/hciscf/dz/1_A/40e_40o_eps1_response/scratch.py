from pyscf import gto, scf, mcscf, grad

mol = gto.M(atom="C 0 0 0; C 0 0 1", basis="ccpvdz", verbose=4)
mf = scf.RHF(mol).run()
mc = mcscf.CASSCF(mf, 8, 8)
mc.kernel()

grad = mc.Gradients().kernel()
print(grad)
