import matplotlib.pyplot as plt
import numpy as np

import uncertainties as unp
from scipy import stats
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
import scipy.constants as const

U_A, delta_y = np.genfromtxt('content/data/data_1.txt', unpack = True)
plt.plot(U_A,delta_y,'bx',label='Lokale Stromänderung')
plt.axvline(x=7, ymin=0, ymax=0.95,color='red')
plt.xlabel(r'$U_A \,/\,\unit{\volt}$')
plt.ylabel(r'$\symup{\Delta}I_A\,/\,\unit{\ampere}$')
plt.legend(loc='best')
plt.grid(which="both")
plt.tight_layout(pad=0, h_pad=1.10, w_pad=1.08)
plt.savefig('build/Differentielle_Energie_20Grad.pdf')
plt.close()

U_A, delta_y = np.genfromtxt('content/data/data_2.txt', unpack = True)
plt.plot(U_A,delta_y,'bx',label='Lokale Stromänderung')
plt.xlabel(r'$U_A \,/\,\unit{\volt}$')
plt.ylabel(r'$\symup{\Delta}I_A\,/\,\unit{\ampere}$')
plt.legend(loc='best')
plt.grid(which="both")
plt.tight_layout(pad=0, h_pad=1.10, w_pad=1.08)
plt.savefig('build/Differentielle_Energie_145Grad.pdf')
plt.close()


E_1 = [5,5,5,5.5]
E_1_mean = np.mean(E_1)
E_1_err = np.std(E_1)
E_1_= ufloat(E_1_mean,E_1_err)
c=3*10**8
h=4.136*10**-15
print(f'Anregungsenergie Reihe 1: {E_1_} eV')
print(f'Emittierte Strahlung 1: {c*h/E_1_*10**9} nm')

E_2 = [4.5,5,5.3,5.1,5.1,5.1,5.9]
E_2_mean = np.mean(E_2)
E_2_err = np.std(E_2)
E_2_= ufloat(E_2_mean,E_2_err)
print(f'Anregungsenergie Reihe 2: {E_2_} eV')
print(f'Emittierte Strahlung 2: {c*h/E_2_*10**9} nm')


def p(T):
    return 5.5*10**7*np.e**(-6876/T)

def w(T):
    return 0.0029/p(T)

def verh(T):
    return 1/w(T)

def korrekturistdoof(T):
    print(f'--------------------')
    print(f'T = {T} °C')
    print(f'p({T}°C) = {p(T+273.15)} mbar')
    print(f'w({T}°C) = {w(T+273.15)} cm')
    print(f'Verhaeltnis: {verh(T+273.15)}')
    print(f'--------------------')
    return

korrekturistdoof(20)
korrekturistdoof(145)
korrekturistdoof(166)
korrekturistdoof(195)
    