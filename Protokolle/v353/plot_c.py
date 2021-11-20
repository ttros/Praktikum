import matplotlib.pyplot as plt
import numpy as np
import uncertainties as unp
from scipy.optimize import curve_fit

f, A, a, b = np.genfromtxt('content/data_bc.txt', unpack = True)
phi = a/b * 2 * np.pi

def f1(f,c):
    return np.arctan(-f*c)

parameters, pcov = curve_fit(f1, f, phi, sigma=None)
print(parameters, np.sqrt(np.diag(pcov)), sep='\n')
plt.plot(np.log(f), phi, 'rx', label='Daten')
plt.plot(np.log(f), f1(f,*parameters), 'b', label='Fit')
#plt.plot(np.log(f), curve_fit(f1,f,A/5.2)))

plt.xlabel(r'log(f) [Hz]')
plt.ylabel(r'phi [rad]')
plt.legend(loc='best')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_c.pdf')

