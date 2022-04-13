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
c_h2o = 1483 #m/s
c_acryl = 2730 #m/s

# %%%%% Schieblehre %%%%%
loch_sch, oben_sch_, unten_sch_ = np.genfromtxt('content/data/schieblehre.txt', unpack = True)
oben_sch = oben_sch_ * 10**-3   #in m
unten_sch = unten_sch_ * 10**-3 #in m

d_sch = hoehe - oben_sch - unten_sch 

print(f'Durchmesser Schieblehre in mm: {d_sch*10**3}')


# %%%%% c-Bestimmung, Kopplungsschicht %%%%%
loch_c, oben_c_ = np.genfromtxt('content/data/c-messung.txt', unpack = True)
oben_sch_kurz = oben_sch
oben_sch_kurz = np.delete(oben_sch_kurz, [7,8,9,10])

oben_c = oben_c_ *10**-6 /2# in micro s, /2 wegen doppelter Weglaenge

plt.plot(oben_sch_kurz,oben_c,'rx',label='Messwerte')

m , b , r ,p ,std =stats.linregress(oben_sch_kurz,oben_c)
M=unp.uarray(m,std)
B=unp.uarray(b,std)

print(f'm: {M}')
print(f'c: {1/M}')
print(f'b: {B}')
print(f'Kopplungsschicht nach Plot in mm: {B*c_h2o*10**3}')
xx = np.linspace(0.01,0.065, 1000)
plt.plot(xx,m*xx+b, 'b', label='Fit')

plt.xlabel(r'$\symup{Tiefe} \,/\,\unit{\metre}')
plt.ylabel(r'$\frac{\symup{\Delta}t}{2}$ \,/\, \unit{\second}')
plt.legend(loc='best')
plt.grid(which="both")

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/c-bestimmung.pdf')
plt.close()

# %%%%% Kopplungsschicht rechnerisch %%%%%
oben_c_theoretisch = oben_sch_kurz / c_acryl # theoretische Laufzeit [s] = tiefe [m] / c [m/s]
#print(f'oben_c_theoretisch: {oben_c_theoretisch*10**6}')
diff_c = oben_c - oben_c_theoretisch
#print(f'diff_c: {diff_c*10**6}')
d_kopplung_array = c_h2o * diff_c
d_kopplung = ufloat(np.mean(d_kopplung_array), np.std(d_kopplung_array))
print(f'Kopplungsschicht nach Rechnung in mm: {d_kopplung*10**3}')
print(f'noms: {noms(d_kopplung)}')

# %%%%% A-Scan %%%%%
loch_a, oben_a_, unten_a_ = np.genfromtxt('content/data/a-scan.txt', unpack = True)
oben_a = oben_a_ * 10**-3 - 2*noms(d_kopplung)  #in m ohne Anpassungsschicht
unten_a = unten_a_ * 10**-3 - 2*noms(d_kopplung)#in m ohne Anpassungsschicht

d_a = hoehe - oben_a - unten_a
d_a_ausgabe = np.around(d_a*10**3, decimals = 1)
print(f'Durchmesser A-Scan      in mm: {d_a*10**3}')

# %%%%% B-Scan %%%%%
loch_b, oben_b_, unten_b_ = np.genfromtxt('content/data/b-scan.txt', unpack = True)
oben_b = oben_b_ * 10**-3  - 2*noms(d_kopplung) #in m ohne Anapssungsschicht
unten_b = unten_b_ * 10**-3 - 2*noms(d_kopplung)#in m ohne Anpassungsschicht

d_b = hoehe - oben_b - unten_b 
d_b_ausgabe = np.around(d_b*10**3, decimals = 2)
print(f'Durchmesser B-Scan      in mm: {d_b*10**3}')
print(f'\n')
print(f'\n')
print(f'\n')
print(f'\n')
print(f'\n')
print(f'Korrigierter A-Scan oben: {np.around(oben_a*10**3, decimals = 1)}')
print(f'Relative Abweichung: {np.around(np.abs(oben_sch-oben_a)/oben_sch*100, decimals = 2)}')
print(f'Korrigierter A-Scan unten: {np.around(unten_a*10**3, decimals = 1)}')
print(f'Relative Abweichung: {np.around(np.abs(unten_sch-unten_a)/unten_sch*100, decimals = 2)}')
print(f'Korrigierter B-Scan oben: {np.around(oben_b*10**3, decimals = 1)}')
print(f'Relative Abweichung: {np.around(np.abs(oben_sch-oben_b)/oben_sch*100, decimals = 2)}')
print(f'Korrigierter B-Scan unten: {np.around(unten_b*10**3, decimals = 1)}')
print(f'Relative Abweichung: {np.around(np.abs(unten_sch-unten_b)/unten_sch*100,decimals = 2)}')
print(f'\n')
print(f'\n')
print(f'\n')
print(f'\n')
print(f'\n')
