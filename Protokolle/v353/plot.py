import matplotlib.pyplot as plt
import numpy as np
import uncertainties as unp
from scipy import stats
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)

t, u = np.genfromtxt('content/data_a.txt', unpack = True)
T = unp.uarray(t,1) #Fehler von t
U = unp.uarray(u,0.1) #Fehler von U

m , b , r ,p ,std =stats.linregress(noms(T),np.log(noms(U/5)))
M=unp.uarray(m,std)
B=unp.uarray(b,std)

plt.plot(t, m*t+b, 'b', label='Fit')
plt.annotate(r'$\symup{ln}(\frac{U}{U_{0}}) =  0,025\cdot t + 1,595$', [-0.8,-0.8])
plt.errorbar(t, np.log(u/5), xerr = stds(T), yerr = stds(U), fmt = "r.", label='Daten')
#plt.plot(t, np.log(U), 'rx', label='Daten')
plt.xlabel(r'$t$ [$\symup{\mu}$s]')
plt.ylabel(r'$\symup{ln}(\frac{U}{U_{0}})$ [V]')
plt.legend(loc='best')

print(f'RC_a = {-1/M} [\mu s]')
print(f'Steigung m: {M}')
print(f'Achsenabschnitt b: {B}')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_a.pdf')
plt.close()

## Beginn Plot B
f, A, a, b = np.genfromtxt('content/data_bc.txt', unpack = True)

w=2*np.pi*f

W_f = unp.uarray(w,0.1) #Fehler von f
A_f = unp.uarray(A,0.1) #Fehler von A

def f1(w,c):
    return 1/(np.sqrt(1+(w**2 * c**2)))

parameters, pcov = curve_fit(f1, w , A/2.8, sigma=None)
RC_B=unp.uarray(parameters,pcov)
#print(f'RC_b = {parameters*10**6} [\mu s]')
print(f'RC_b = {RC_B*10**6} [\mu s]')
#plt.plot(np.log(f), A/2.8, 'rx', label='Daten')
plt.errorbar(np.log(w), A/2.8 , xerr = stds(W_f), yerr = stds(A_f), fmt = 'r.', label='Daten')
xx = np.linspace((np.e)**5, (np.e)**14, 10000)
plt.plot(np.log(xx), f1(xx,*parameters), 'b', label='Fit')

plt.xlabel(r'$\symup{ln}(\omega)$ [Hz]')
plt.ylabel(r'$A/U_{0}$')
plt.legend(loc='best')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_b.pdf')
plt.close()

## Beginn Plot C

a_f = unp.uarray(a,0.1) #Fehler von a
b_f = unp.uarray(b,0.1) #Fehler von b

phi = a_f/b_f * 2 * np.pi

def f2(w,c):
    return np.arctan(-w*c)

parameters, pcov = curve_fit(f2, w, noms(phi), sigma=None)
RC_C=unp.uarray(parameters,pcov)
print(f'RC_c = {-RC_C*10**6} [\mu s]')
#plt.plot(np.log(f), phi, 'rx', label='Daten')
plt.errorbar(np.log(w), noms(phi) , xerr = stds(W_f), yerr = stds(phi), fmt = 'r.', label='Daten')
plt.plot(np.log(xx), f2(xx,*parameters), 'b', label='Fit')

plt.xlabel(r'$\symup{log}(\omega)$ [Hz]')
plt.ylabel(r'$\varphi$ [rad]')
plt.legend(loc='best')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_c.pdf')
plt.close()

## Beginn Plot D
phi2 = np.linspace(0, (np.pi)/2, 50)

plt.polar(noms(phi), A, 'rx')
#plt.polar(phi2, np.sin(phi2)/(np.tan(phi2))*2.8, 'b')
plt.polar(phi2, np.cos(phi2)*2.8, 'b')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_d.pdf')
plt.close()
