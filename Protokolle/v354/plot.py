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

#R=480????????????????????????????
L=10.11*10**-3
C=2.093*10**-9

f_kHz, U_0, U, a, b = np.genfromtxt('content/data/data_c_d.txt', unpack = True)

f = f_kHz*1000
U_RCL = U/U_0

#def Theorie_c(f, R, L, C):
#    return (1/np.sqrt((1-L*C*f**2)**2+f**2*R**2*C**2))

plt.plot(f_kHz, U_RCL, 'r', label='Messwertekurve')
#plt.plot(f_kHz, Theorie_c(2*np.pi*1000*f_kHz, R, L, C), label = "Theoriekurve")
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

plt.plot(f_kHz, U_RCL, 'r', label='Messwertekurve')

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

plt.plot(f_kHz, delta_phi, 'r', label='Messwertekurve')
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

plt.plot(f_kHz, delta_phi, 'r', label='Messwertekurve')

plt.xlabel(r'$f/ \unit{\kilo\hertz}$')
plt.ylabel(r'$\symup{\Delta}\varphi$')
plt.legend(loc='best')
plt.grid(which="both")
plt.ylim(0.75 * np.pi/2, 1.25 * np.pi/2)
plt.yticks([0, np.pi / 4, np.pi / 2, 3 * np.pi/4, np.pi],
           [r"$0$", r"$\frac{1}{4}\pi$", r"$\frac{1}{2}\pi$", r"$\frac{3}{4}\pi$", r"$\pi$"])

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_d_2.pdf')
plt.close()