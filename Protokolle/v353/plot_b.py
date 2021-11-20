import matplotlib.pyplot as plt
import numpy as np
import uncertainties as unp
from scipy.optimize import curve_fit

f, A, a, b = np.genfromtxt('content/data_bc.txt', unpack = True)


def f1(f,c):
    return 1/(np.sqrt(1+(f**2 * c**2)))

parameters, pcov = curve_fit(f1, f, A/2.8, sigma=None)
print(parameters, np.sqrt(np.diag(pcov)), sep='\n')
plt.plot(np.log(f), A/2.8, 'rx', label='Daten')
plt.plot(np.log(f), f1(f,*parameters), 'b', label='Fit')

plt.xlabel(r'log(f) [Hz]')
plt.ylabel(r'A/U_0')
plt.legend(loc='best')

# # in matplotlibrc leider (noch) nicht möglich
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_b.pdf')

