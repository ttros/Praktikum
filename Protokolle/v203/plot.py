import matplotlib.pyplot as plt
import numpy as np

import uncertainties as unp
from scipy import stats
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)

# %%%%% Daten 1 in SI-Einheiten %%%%% #
T_messung, p_messung = np.genfromtxt('content/data/data_1.txt', unpack = True)


T = T_messung + 273.15      # T in K
p = p_messung * 10**(-1)       # p in kPa

np.savetxt('data_1_si.txt', np.column_stack([T, p]),fmt=['%d', '%.1f'], header="T / K, p / kPa")
# %%%%% Daten 1 in SI-Einheiten %%%%% #

############### Plot 1 ################
plt.plot(T,p,'rx', label='Messwerte')

def f(T,A,B):
    return A*(np.e)**(-B*1/T)

parameters, pcov = curve_fit(f, T, p, sigma=None)
A = parameters[0]
B = parameters[1]
print(f'\n\nA: {A}')
print(f'B: {B}')
print(f'L = {B*8.314}')
print(f'\n\n')

xx = np.linspace(20, 100, 10000)
plt.plot(T, f(T,*parameters), 'b', label = 'Fit')

#plt.yscale('log')

plt.xlabel(r'$T \,/\,\unit{\kelvin}$')
plt.ylabel(r'$p \,/\, \unit{\kilo\pascal}$')
plt.legend(loc='best')
plt.grid(which="both")

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_1.pdf')
