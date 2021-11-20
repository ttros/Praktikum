import numpy as np
s, U, a, b = np.genfromtxt('data_cd.txt', unpack = True)

A_U = U/5.2
phi = a/b * 2 * np.pi
np.savetxt('A_U.txt', np.column_stack([A_U]), fmt='%.3f', header='A/U_0')

