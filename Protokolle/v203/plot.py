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

# %%%%%%%%%%%%%%%%%%%%%% TEIL 1 %%%%%%%%%%%%%%%%%%%%%% #
####### Daten 1 in SI-Einheiten #######
T_messung, p_messung = np.genfromtxt('content/data/data_1.txt', unpack = True)


T = T_messung + 273.15      # T in K
p = p_messung * 10**(-1)    # p in kPa
p_0 = 101.7                 # p_0 in kPa 

np.savetxt('content/data/data_1_si.txt', np.column_stack([T, p]),fmt=['%d', '%.1f'], header="T / K, p / kPa")
####### Daten 1 in SI-Einheiten #######

############### Plot 1 ################
plt.plot(1/T,np.log(p/p_0),'rx', label='Messwerte')

m , b , r ,p ,std =stats.linregress(1/T,np.log(p/p_0))
M=unp.uarray(m,std)
B=unp.uarray(b,std)

plt.plot(1/T, m*1/T +b, 'b', label = 'Fit')

plt.xlabel(r'$\frac{1}{T} \,/\,\frac{1}{\unit{\kelvin}}$')
plt.ylabel(r'$\ln\left(\frac{p}{p_0}\right)$')
plt.legend(loc='best')
plt.grid(which="both")

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_1.pdf')
plt.close()

############### Plot 1 ################
############ Berechnung L #############
L = -M*8.314
L_i = L - 3101
L_i_m_temp = L/(6.02*10**23) # Einheit J
L_i_m = L_i_m_temp * 6.242*10**18 # Einheit eV

print(f'\n\n')
print(f'M: {M} K')
print(f'B: {B}')
print(f'L: {L} J/(mol)')
print(f'L_i: {L_i} J/(mol)')
print(f'L_i_m: {L_i_m} eV')
print(f'\n\n')
############ Berechnung L #############
# %%%%%%%%%%%%%%%%%%%%%% TEIL 1 %%%%%%%%%%%%%%%%%%%%%% #
# %%%%%%%%%%%%%%%%%%%%%% TEIL 2 %%%%%%%%%%%%%%%%%%%%%% #
####### Daten 2 in SI-Einheiten #######
p_messung_2, T_messung_2 = np.genfromtxt('content/data/data_2.txt', unpack = True)

p_2 = p_messung_2 * 10**(2)    # p in kPa
T_2 = T_messung_2 + 273.15      # T in K

np.savetxt('content/data/data_2_si.txt', np.column_stack([p_2, T_2]),fmt=['%d', '%d'], header="p / kPa, T / K") 
####### Daten 2 in SI-Einheiten #######
############### Plot 2 ################
plt.plot(T_2,p_2,'rx', label='Messwerte')

def f(T,a,b,c,d):
    return a*(T**3) + b*(T**2) + c*T + d

parameters, pcov = curve_fit(f, T_2 , p_2, sigma=None)
std = np.sqrt(np.diag(pcov))

xx = np.linspace(388,472, 1000)
plt.plot(xx, f(xx,*parameters), 'b', label='Fit')

plt.xlabel(r'$T \,/\,\unit{\kelvin}$')
plt.ylabel(r'$p \,/\, \unit{\kilo\pascal}$')
plt.legend(loc='best')
plt.grid(which="both")

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_2.pdf')
plt.close()
############### Plot 2 ################
############## Ausgabe 2 ##############
A = ufloat(parameters[0], std[0])
B = ufloat(parameters[1], std[1])
C = ufloat(parameters[2], std[2])
D = ufloat(parameters[3], std[3])
print(f'Parameter Fit 2. Teil')
print(f'A: {A}')
print(f'B: {B}')
print(f'C: {C}')
print(f'D: {D}')
print(f'\n\n')
############## Ausgabe 2 ##############
############## Plot 3,4 ###############
R=const.gas_constant
C=0.9
def L_plus(x,a,b,c,d):
    return ((((R*x)/2)+np.sqrt(((R**2*x**2)/2)-C*(a*x**3+b*x**2+c*x+d)))*((3*a*x**3+2*b*x**2+c*x)/(a*x**3+b*x**2+c*x+d)))
def L_minus(x,a,b,c,d):
    return ((((R*x)/2)-np.sqrt(((R**2*x**2)/2)-C*(a*x**3+b*x**2+c*x+d)))*((3*a*x**3+2*b*x**2+c*x)/(a*x**3+b*x**2+c*x+d)))

x = np.linspace(390, 470, 1000)

plt.plot(x, L_plus(x,*parameters*10**3), 'r')

plt.xlabel(r'$T \,/\,\unit{\kelvin}$')
plt.ylabel(r'$L_{+}(T) \,/\, \unit{\joule\per\mol}$')
plt.grid(which="both")

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_L_plus.pdf')
plt.close()



plt.plot(xx, L_minus(x,*parameters*10**3), 'r')

plt.xlabel(r'$T \,/\,\unit{\kelvin}$')
plt.ylabel(r'$L_{-}(T) \,/\, \unit{\joule\per\mol}$')
plt.grid(which="both")

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_L_minus.pdf')
plt.close()
############## Plot 3,4 ###############
# %%%%%%%%%%%%%%%%%%%%%% TEIL 2 %%%%%%%%%%%%%%%%%%%%%% #