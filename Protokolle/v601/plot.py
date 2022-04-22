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