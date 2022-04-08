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
l_, t_ = np.genfromtxt('content/data/c-Bestimmung_IE.txt', unpack = True)
l = l_ * 10**(-3)
t = t_/2 * 10**(-6)
plt.plot(l[-1],t[-1], color='dimgray', marker='x')
l = np.delete(l, -1)
t = np.delete(t, -1)
print(t[-2])
plt.plot(l,t,'rx',label='Messwerte')

m , b , r ,p ,std =stats.linregress(l,t)
M=unp.uarray(m,std)
B=unp.uarray(b,std)

print(f'm: {M}')
print(f'c: {1/M}')
print(f'b: {B}')

xx = np.linspace(38*10**-3,162*10**-3, 1000)
plt.plot(xx,m*xx+b, 'b', label='Fit')

plt.xlabel(r'$l \,/\,\unit{\metre}')
plt.ylabel(r'$\symup{\Delta}t$ \,/\, \unit{\second}')
plt.legend(loc='best')
plt.grid(which="both")

plt.tight_layout(pad=0, h_pad=1.10, w_pad=1.08)
plt.savefig('build/c-Bestimmung_IE.pdf')
plt.close()

# %%%%% c-Bestimmung D-Verfahren %%%%%
l_2, t_2 = np.genfromtxt('content/data/c-Bestimmung_D.txt', unpack = True)
l2 = l_2 * 10**(-3)
t2 = t_2 * 10**(-6)

plt.plot(l2,t2,'rx',label='Messwerte')

m , b , r ,p ,std =stats.linregress(l,t)
M=unp.uarray(m,std)
B=unp.uarray(b,std)

print(f'm2: {M}')
print(f'c2: {1/M}')
print(f'b2: {B}')

xx = np.linspace(38*10**-3,122*10**-3, 1000)
plt.plot(xx,m*xx+b, 'b', label='Fit')

plt.xlabel(r'$l \,/\,\unit{\metre}')
plt.ylabel(r'$\symup{\Delta}t$ \,/\, \unit{\second}')
plt.legend(loc='best')
plt.grid(which="both")

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/c-Bestimmung_D.pdf')
plt.close()

# %%%%% Daempfungsbestimmung IE-Verfahren %%%%%
l_3, A_out, A_in, output = np.genfromtxt('content/data/Daempfung_IE.txt', unpack = True)

A = A_in/A_out
l3= 2*l_3 * 10**-3
print(f'A: {A}')
plt.plot(l3[-1],np.log(A[-1]), color='dimgray', marker='x')
plt.plot(l3[-2],np.log(A[-2]), color='dimgray', marker='x')
l3 = np.delete(l3, [5,6])
A = np.delete(A, [5,6])
print(l3)
plt.plot(l3,np.log(A),'rx',label='Messwerte') #Achsenbeschriftung anpassen!

m , b , r ,p ,std =stats.linregress(l3,np.log(A))
M=unp.uarray(m,std)
B=unp.uarray(b,std)

print(f'm2: {M}')
print(f'c2: {1/M}')
print(f'b2: {B}')

xx = np.linspace(2*38*10**-3,2*162*10**-3, 1000)
plt.plot(xx,m*xx+b, 'b', label='Fit')

plt.xlabel(r'$2 \cdot l \,/\,\unit{\metre}')
plt.ylabel(r'$\ln(\frac{A_{\symup{in}}}{A_{\symup{out}}})$')
plt.legend(loc='best')
plt.grid(which="both")

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Daempfungsbestimmung.pdf')
plt.close()

# %%%%% Abmessungen Auge %%%%%
a, t_4 = np.genfromtxt('content/data/Auge.txt', unpack = True)
t4 = t_4 / 2 * 10**-6
print(f't4: {t4}')
c_1 = 1483
c_2 = 2500
c_3 = 1410

def d(t,c,c_old):
    return (c*t*10**2 - c_old)  #output in cm

d1 = d(t4[0],c_1,0)
d2 = d(t4[1],c_1,d1)
d3 = d(t4[2],c_2,d1+d2)
d4 = d(t4[3],c_3,d1+d2+d3)
print(f'd1: {d1}')
print(f'd2: {d2}')
print(f'd3: {d3}')
print(f'd4: {d4}')