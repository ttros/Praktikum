import matplotlib.pyplot as plt
import numpy as np

import uncertainties as unp
from scipy import stats
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)

phi, U_normal, U_rausch = np.genfromtxt('content/Daten/data.txt', unpack = True)

phi_rad = phi*2*np.pi/360 

def f(phi_rad,A,B):
    return np.sqrt((A*np.cos(phi_rad))**2)+B

# Amplitude normal

parameters, pcov = curve_fit(f, phi_rad , U_normal, sigma=None)
std = np.sqrt(np.diag(pcov))
print(std)
A= ufloat(parameters[0], std[0])
B= ufloat(parameters[1], std[1])
print(f'Amplitude normal = {A} [\V]')
plt.plot(phi_rad, U_normal, 'rx', label='Daten')
xx = np.linspace(0, 2*np.pi, 10000)
plt.plot(xx, f(xx,*parameters), 'b', label='Fit')

plt.xlabel(r'$\phi / \unit{\degree}$')
plt.ylabel(r'$U / \unit{\volt}$')
plt.legend(loc='best')
plt.grid()

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_normal.pdf')
plt.close()

# LED Messung

r, U_led = np.genfromtxt('content/Daten/data_led.txt', unpack = True)

def f_led(r,C,M):
    return C/(r**2)+M

parameters, pcov = curve_fit(f_led, r , U_led, sigma=None)
std = np.sqrt(np.diag(pcov))
C= ufloat(parameters[0], std[0])
print(f'C bei 1/r^2 = {C} [\V]')
plt.plot(r, U_led, 'rx', label='Daten')
xx = np.linspace(6, 60, 10000)
plt.plot(xx, f_led(xx,*parameters), 'b', label='Fit')

plt.xlabel(r'$r / \unit{\centi\meter}$')
plt.ylabel(r'$U / \unit{\volt}$')
plt.legend(loc='best')
plt.grid()

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_led.pdf')
plt.close()