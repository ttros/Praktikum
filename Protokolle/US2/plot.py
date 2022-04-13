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

hoehe = 80.5 * 10**-3 #in m

# %%%%% Schieblehre %%%%%
loch_sch, oben_sch_, unten_sch_ = np.genfromtxt('content/data/schieblehre.txt', unpack = True)
oben_sch = oben_sch_ * 10**-3   #in m
unten_sch = unten_sch_ * 10**-3 #in m

d_sch = hoehe - oben_sch - unten_sch 

print(f'Durchmesser Schieblehre in mm: {d_sch*10**3}')


# %%%%% c-Bestimmung %%%%%
loch_c, oben_c = np.genfromtxt('content/data/c-messung.txt', unpack = True)
oben_sch_kurz = oben_sch
oben_sch_kurz = np.delete(oben_sch_kurz, [7,10])

print(f'DA: {oben_sch_kurz}')

# plt.plot(l2,t2,'rx',label='Messwerte')

# m , b , r ,p ,std =stats.linregress(l,t)
# M=unp.uarray(m,std)
# B=unp.uarray(b,std)

# xx = np.linspace(38*10**-3,122*10**-3, 1000)
# plt.plot(xx,m*xx+b, 'b', label='Fit')

# plt.xlabel(r'$l \,/\,\unit{\metre}')
# plt.ylabel(r'$\symup{\Delta}t$ \,/\, \unit{\second}')
# plt.legend(loc='best')
# plt.grid(which="both")

# plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
# plt.savefig('build/c-bestimmung.pdf')
# plt.close()

# %%%%% A-Scan %%%%%
loch_a, oben_a_, unten_a_ = np.genfromtxt('content/data/a-scan.txt', unpack = True)
oben_a = oben_a_ * 10**-3   #in m
unten_a = unten_a_ * 10**-3 #in m

d_a = hoehe - oben_a - unten_a +1.19*2*10**-3

print(f'Durchmesser A-Scan      in mm: {d_a*10**3}')

# %%%%% B-Scan %%%%%
loch_b, oben_b_, unten_b_ = np.genfromtxt('content/data/b-scan.txt', unpack = True)
oben_b = oben_b_ * 10**-3   #in m
unten_b = unten_b_ * 10**-3 #in m

d_b = hoehe - oben_b - unten_b +1.19*2*10**-3

print(f'Durchmesser B-Scan      in mm: {d_b*10**3}')

###############
x = np.linspace(0, 10, 1000)
y = x ** np.sin(x)

plt.subplot(1, 2, 1)
plt.plot(x, y, label='Kurve')
plt.xlabel(r'$\alpha \mathbin{/} \unit{\ohm}$')
plt.ylabel(r'$y \mathbin{/} \unit{\micro\joule}$')
plt.legend(loc='best')

plt.subplot(1, 2, 2)
plt.plot(x, y, label='Kurve')
plt.xlabel(r'$\alpha \mathbin{/} \unit{\ohm}$')
plt.ylabel(r'$y \mathbin{/} \unit{\micro\joule}$')
plt.legend(loc='best')

# in matplotlibrc leider (noch) nicht möglich
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot.pdf')
