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

# %%%%% c-Bestimmung I-E-Verfahren %%%%%
l, t_ = np.genfromtxt('content/data/c-Bestimmung_IE.txt', unpack = True)
t = t_/2
print(t[-2])
plt.plot(l,t,'rx',label='Messwerte')

plt.xlabel(r'$l \,/\,\unit{\milli\metre}')
plt.ylabel(r'$\symup{\Delta}t$ \,/\, \unit{\micro\second}')
plt.legend(loc='best')
plt.grid(which="both")

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/c-Bestimmung_IE.pdf')
plt.close()

# %%%%% c-Bestimmung D-Verfahren %%%%%
l, t = np.genfromtxt('content/data/c-Bestimmung_D.txt', unpack = True)
plt.plot(l,t,'rx',label='Messwerte')

plt.xlabel(r'$l \,/\,\unit{\milli\metre}')
plt.ylabel(r'$\symup{\Delta}t$ \,/\, \unit{\micro\second}')
plt.legend(loc='best')
plt.grid(which="both")

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/c-Bestimmung_D.pdf')
plt.close()

# %%%%% Daempfungsbestimmung IE-Verfahren %%%%%
l, A_out, A_in, output = np.genfromtxt('content/data/Daempfung_IE.txt', unpack = True)

A = A_in/A_out

plt.plot(l,np.log(A),'rx',label='Messwerte') #Achsenbeschriftung anpassen!


plt.xlabel(r'$l \,/\,\unit{\milli\metre}')
plt.ylabel(r'$\symup{\Delta}t$ \,/\, \unit{\micro\second}')
plt.legend(loc='best')
plt.grid(which="both")

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Daempfungsbestimmung.pdf')
plt.close()