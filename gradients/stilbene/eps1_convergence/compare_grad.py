import numpy as np

g_103 = np.load("_data/vhciscf_1.0e-03.npy")
g_504 = np.load("_data/vhciscf_5.0e-04.npy")
g_104 = np.load("_data/vhciscf_1.0e-04.npy")
g_105 = np.load("_data/vhciscf_1.0e-05.npy")

print(np.linalg.norm(g_103-g_105))
print(np.std((g_103-g_105).flatten()))
print(np.std((g_104-g_105).flatten()))