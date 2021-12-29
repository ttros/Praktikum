import matplotlib.pyplot as plt
import numpy as np

import uncertainties as unp
from scipy import stats
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)

######################## Plot A #############################
A = np.genfromtxt('content/data/data_a.txt', unpack = True)
t = 27.27 # 27.27in micro sec

x = np.linspace(0,13,14)
xx = np.linspace(0,14, 1000)
plt.subplot(1, 2, 1)
plt.plot(t*x, A, 'rx', label='Messwerte')

plt.xlabel(r'$t\,/\,\unit{\micro\second}$')
plt.ylabel(r'$U \,/\, \unit{\volt}$')
plt.legend(loc='best')
plt.grid(which="both")

plt.subplot(1, 2, 2)
m , b , r ,p ,std =stats.linregress(t*x,np.log(noms(A/6)))
M=unp.uarray(m,std)
B=unp.uarray(b,std)
print(f'\n')
print(f'Steigung m: {M}')
print(f'Achsenabschnitt b: {B}')
plt.plot(t*x,np.log(A/6), 'rx', label='Messwerte')
plt.plot(t*xx,m*t*xx+b, 'b', label='Fit')

plt.xlabel(r'$t\,/\,\unit{\micro\second}$')
plt.ylabel(r'$\symup{ln}(\frac{U}{U_{0}})$')
plt.legend(loc='best')
plt.grid(which="both")

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_a.pdf')
plt.close()
######################## Plot C #############################

R=559.5 #Innenwiderstand 50
L=10.11*10**-3 #Gerätedaten
C=2.093*10**-9 #Gerätedaten
f_kHz, U_0, U, a, b = np.genfromtxt('content/data/data_c_d.txt', unpack = True)

f = f_kHz*1000
U_RCL = U/U_0

def Theorie_c(w, R, L, C):
    return (1/np.sqrt((1-L*C*w**2)**2+w**2*R**2*C**2))

x = np.linspace(10, 47, 1000)

plt.plot(f_kHz, U_RCL, 'rx', label='Messwerte')
plt.plot(x, Theorie_c(2*np.pi*1000*x, R, L, C), label = "Theoriekurve")
plt.plot(33, 14/3.6 , 'go', label = "Maximum der Messwerte")
plt.xscale('log')

plt.xlabel(r'$f/ \unit{\kilo\hertz}$')
plt.ylabel(r'$\frac{U_0}{U}$')
plt.legend(loc='best')
plt.grid(which="both")

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_c.pdf')
plt.close()

## plot 2 ##

def Breite(x):
    return (3.9/np.sqrt(2))*x/x

x_b = np.linspace(28.9, 37.7, 2)

plt.plot(f_kHz, U_RCL, 'r', label='Messwertekurve')
plt.plot(x_b, Breite(x_b),'c--', label = r'Breite der Messwertekurve: $\symup{\Delta} f = 8.8\,\unit{\kilo\hertz}$')
plt.plot(33, 14/3.6 , 'go', label = "Maximum der Messwerte")

plt.xlabel(r'$f/ \unit{\kilo\hertz}$')
plt.ylabel(r'$\frac{U_0}{U}$')
plt.legend(loc='best')
plt.grid(which="both")
plt.ylim(0, 5)
plt.xlim(25, 40)

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_c_2.pdf')
plt.close()

######################## Plot D #############################

delta_phi = (a/b) * 2 * np.pi

def Theorie_d(w, R, L, C):
    return np.pi/2+np.arctan((1-L*C*(w)**2)/(-w*R*C))

plt.plot(f_kHz, delta_phi, 'rx', label='Messwerte')
plt.plot(x, Theorie_d(2*np.pi*1000*x, R, L, C), label='Theoriekurve')

plt.xscale('log')

plt.xlabel(r'$f/ \unit{\kilo\hertz}$')
plt.ylabel(r'$\symup{\Delta}\varphi$')
plt.legend(loc='best')
plt.grid(which="both")
plt.yticks([0, np.pi / 4, np.pi / 2, 3 * np.pi/4, np.pi],
           [r"$0$", r"$\frac{1}{4}\pi$", r"$\frac{1}{2}\pi$", r"$\frac{3}{4}\pi$", r"$\pi$"])

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_d.pdf')
plt.close()

## plot 2 ##

plt.plot(f_kHz, delta_phi, 'r', label='Messwertekurve')
plt.axvline(30.4, color='tab:orange', linestyle=':', label=r'$f_1 = 30.4\,\unit{\kilo\hertz}$ und $f_2 = 39.4\,\unit{\kilo\hertz}$')
plt.axvline(39.4, color='tab:orange', linestyle=':')
plt.axvline(34.3, color='g', linestyle=':', label=r'$f_{\symup{res}} = 34.3\,\unit{\kilo\hertz}$')

plt.xlabel(r'$f/ \unit{\kilo\hertz}$')
plt.ylabel(r'$\symup{\Delta}\varphi$')
plt.legend(loc='best')
plt.grid(which="both")
plt.xlim(28, 42)
plt.yticks([0, np.pi / 4, np.pi / 2, 3 * np.pi/4, np.pi],
           [r"$0$", r"$\frac{1}{4}\pi$", r"$\frac{1}{2}\pi$", r"$\frac{3}{4}\pi$", r"$\pi$"])

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_d_2.pdf')
plt.close()

####### Theorie Berechnungen #######

R1_f = ufloat(48.1,0.1)
R2_f = ufloat(509.5,0.5)
L_f = ufloat(10.11*(10**(-3)),0.03*(10**(-3)))
C_f = ufloat(2.093*(10**(-9)),0.003*(10**(-9)))

print(f'\n')

## a ##
R_eff = R1_f + 50
print(f'R effektiv Theorie: {R_eff} Ohm')
print(f'R effektiv Experiement: {-2*(M*10**6)*L} Ohm')

## b ##
R_ap = 2*unp.sqrt(L_f/C_f)
print(f'R aperiodisch: {R_ap/1000} kOhm')

## c ##
delta_omega = (R2_f+50)/L
delta_f = delta_omega/(2*np.pi)
omega_0 = unp.sqrt(1/(L*C))
q = omega_0/delta_omega
print(f'Breite der Resonanzkurve: {delta_f/1000} kHz')
print(f'Güte bzw. Resonanzüberhöhung: {q} V')

## d ##
omega_res = unp.sqrt(1/(L*C) - ((R2_f+50)**2)/(2*(L**2)))
omega_1 = (R2_f+50)/(2*L) + unp.sqrt( (R2_f+50)**2/(4*(L**2)) + 1/(L*C) )
omega_2 = -(R2_f+50)/(2*L) + unp.sqrt( (R2_f+50)**2/(4*(L**2)) + 1/(L*C) )
print(f'Resonanzfrequenz: {omega_res/(1000*2*np.pi)} kHz')
print(f'Frequenz 1: {omega_1/(1000*2*np.pi)} kHz')
print(f'Frequenz 2: {omega_2/(1000*2*np.pi)} kHz')

######
print(f'\n')

##### Ausgabe #####

#Steigung m: -0.0106+/-0.0007
#Achsenabschnitt b: 0.3013+/-0.0007
#
#
#R effektiv Theorie: 98.10+/-0.10 Ohm
#R effektiv Experiement: 214+/-14 Ohm
#R aperiodisch: 4.396+/-0.007 kOhm
#Breite der Resonanzkurve: 8.808+/-0.008 kHz
#Güte bzw. Resonanzüberhöhung: 3.9282+/-0.0035 V
#Resonanzfrequenz: 34.0335+/-0.0010 kHz
#Frequenz 1: 39.282+/-0.004 kHz
#Frequenz 2: 30.4739+/-0.0034 kHz