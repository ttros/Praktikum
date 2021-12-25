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

R=480
L=10.11*10**-3
C=2.093*10**-9

f_kHz, U_0, U, a, b = np.genfromtxt('content/data/data_c_d.txt', unpack = True)

f = f_kHz*1000
U_RCL = U/U_0

def Theorie_c(f, R, L, C):
    return (1/np.sqrt((1-L*C*f**2)**2+f**2*R**2*C**2))

plt.plot(f_kHz, U_RCL, 'rx', label='Daten')
plt.plot(f_kHz, Theorie_c(2*np.pi*1000*f_kHz, R, L, C), label = "Theoriekurve")

plt.xscale('log')

plt.xlabel(r'$f/ \unit{\kilo\hertz}$')
plt.ylabel(r'$\frac{U_0}{U}$')
plt.legend(loc='best')
plt.grid(which="both")

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_c.pdf')
plt.close()


######################## Plot D #############################

delta_phi = (a/b) * 2 * np.pi

plt.plot(f_kHz, delta_phi, 'rx', label='Daten')
plt.xscale('log')

plt.xlabel(r'$f/ \unit{\kilo\hertz}$')
plt.ylabel(r'$\symup{\Delta}\varphi$')
plt.legend(loc='best')
plt.grid(which="both")

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_d.pdf')
plt.close()

#def f(phi_rad,A,B):
#    return np.sqrt((A*np.cos(phi_rad))**2)+B
#
#parameters, pcov = curve_fit(f, phi_rad , U_normal, sigma=None)
#std = np.sqrt(np.diag(pcov))
#print(std)
#A= ufloat(parameters[0], std[0])
#B= ufloat(parameters[1], std[1])
#print(f'Amplitude normal = {A} [\V]')
#plt.plot(phi_rad, U_normal, 'rx', label='Daten')
#xx = np.linspace(0, 2*np.pi, 10000)
#plt.plot(xx, (2*48/np.pi)*np.cos(xx), 'b', label='Theorie')
#plt.xticks([0, np.pi / 2, np.pi, 3 / 2 * np.pi, 2 * np.pi],
#          [r"$0$", r"$\frac{1}{2}\pi$", r"$\pi$", r"$\frac{3}{2}\pi$", r"$2\pi$"])
#
#plt.xlabel(r'$\symup{\Delta}\varphi')
#plt.ylabel(r'$U / \unit{\volt}$')
#plt.legend(loc='best')
#plt.grid()
#
#plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
#plt.savefig('build/plot_normal.pdf')
#plt.close()