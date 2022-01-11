import matplotlib.pyplot as plt
import numpy as np

import uncertainties as unp
from scipy import stats
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)

############### Plot 1 ###############
T, p = np.genfromtxt('content/data/data_1.txt', unpack = True)

plt.plot(T,p,'rx', label='Messwerte')

def f(A,B,T):
    return A*(np.e)**(B*T)

parameters, pcov = curve_fit(f, T , p, sigma=None)
A = parameters[0]
B = parameters[1]
print(f'A: {A}')
print(f'B: {B}')

xx = np.linspace(20, 100, 10000)
#plt.plot(xx, f(A,B,xx), 'b', label = 'Fit')

#plt.yscale('log')

plt.xlabel(r'$T\,/\,\unit{\celsius}$')
plt.ylabel(r'$p \,/\, \unit{\milli\bar}$')
plt.legend(loc='best')
plt.grid(which="both")

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_1.pdf')
