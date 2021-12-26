import matplotlib.pyplot as plt
import numpy as np

import uncertainties as unp
from scipy import stats
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)

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

#Plot zwai

plt.plot(f_kHz, U_RCL, 'rx', label='Messwerte')
plt.plot(x, Theorie_c(2*np.pi*1000*x, R, L, C), label = "Theoriekurve")

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
    return  np.arctan((-w*R*C)/(1-L*C*(w**2)))

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

#Plot zwai

plt.plot(f_kHz, delta_phi, 'rx', label='Messwerte')

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