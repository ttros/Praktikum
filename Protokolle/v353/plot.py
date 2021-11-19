import matplotlib.pyplot as plt
import numpy as np
import uncertainties as unp
from scipy.optimize import curve_fit

t, U = np.genfromtxt('content/data_a.txt', unpack = True)

m,b = np.polyfit(t, np.log(U), 1)
plt.plot(t, m*t+b, 'b', label='Fit')
plt.annotate(f'$ln(U) =  {m} \cdot f + {b}$', [0.8,0.2])
#plt.errorbar(t, np.log(U), xerr = 1, yerr = np.e**U * 0.1, fmt = 'rx') #Fehler y anpassen!!! Gaussssss
plt.plot(t, np.log(U), 'rx', label='Daten')
plt.xlabel(r'$f$ [Hz]')
plt.ylabel(r'$ln(U)$ [V]')
plt.legend(loc='best')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_a.pdf')
plt.close()

## Beginn Plot B
f, A, a, b = np.genfromtxt('content/data_bc.txt', unpack = True)
phi = a/b * 2 * np.pi

def f1(f,c):
    return 1/(np.sqrt(1+(f**2 * c**2)))

parameters, pcov = curve_fit(f1, f, A/2.8, sigma=None)
print(parameters, np.sqrt(np.diag(pcov)), sep='\n')
plt.plot(np.log(f), A/2.8, 'rx', label='Daten')
plt.plot(np.log(f), f1(f,*parameters), 'b', label='Fit')

plt.xlabel(r'$f$ [Hz]')
plt.ylabel(r'$A/U_{0}$')
plt.legend(loc='best')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_b.pdf')
plt.close()

## Beginn Plot C
def f2(f,c):
    return np.arctan(-f*c)

parameters, pcov = curve_fit(f2, f, phi, sigma=None)
print(parameters, np.sqrt(np.diag(pcov)), sep='\n')
plt.plot(np.log(f), phi, 'rx', label='Daten')
plt.plot(np.log(f), f2(f,*parameters), 'b', label='Fit')

plt.xlabel(r'$log(f)$ [Hz]')
plt.ylabel(r'$\varphi$ [rad]')
plt.legend(loc='best')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_c.pdf')
plt.close()

## Beginn Plot D
phi2 = np.linspace(0.1, (np.pi)/2, 50)

plt.polar(phi, A, 'rx')
plt.polar(phi2, np.sin(phi2)/(np.tan(phi2))*2.8, 'b')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_d.pdf')
plt.close() 

