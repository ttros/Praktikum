import matplotlib.pyplot as plt
import numpy as np
import uncertainties as unp
from scipy import stats
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)

T_, t_O_1_, t_O_2_, t_U_1_, t_U_2_ = np.genfromtxt('content/data_T.txt', unpack = True)
# Messunsicherheiten
T = unp.uarray(T_,1)
t_O_1 = unp.uarray(t_O_1_,0.5)
t_O_2 = unp.uarray(t_O_2_,0.1)
t_U_1 = unp.uarray(t_U_1_,0.1)
t_U_2 = unp.uarray(t_U_2_,0.1)

# m , b , r ,p ,std =stats.linregress(noms(T),noms(t_O_1))
# M=unp.uarray(m,std)
# B=unp.uarray(b,std)
#plt.plot(t, m*t+b, 'b', label='Fit')
#plt.annotate(f'$ln(U) =  {m} \cdot t + {b}$', [0,0.15])
plt.errorbar(1/noms(T), np.log(1/noms(t_O_1)), xerr = stds(1/T), yerr = stds(t_O_1), fmt = "r.", label='Daten')
#plt.plot(t, np.log(U), 'rx', label='Daten')
#plt.xlabel(r'$t$ [$\symup{\mu}$s]')
#plt.ylabel(r'$ln(U)$ [V]')
#plt.legend(loc='best')
# print(f'RC_a = {-1/M} [\mu s]')
# print(f'Steigung m: {M}')
# print(f'Achsenabschnitt b: {B}')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('plot_a.pdf')
plt.close()

# ## Beginn Plot B
# f, A, a, b = np.genfromtxt('content/data_bc.txt', unpack = True)

# w=2*np.pi*f

# W_f = unp.uarray(w,0.1) #Fehler von f
# A_f = unp.uarray(A,0.1) #Fehler von A

# def f1(w,c):
#     return 1/(np.sqrt(1+(w**2 * c**2)))

# parameters, pcov = curve_fit(f1, w , A/2.8, sigma=None)
# RC_B=unp.uarray(parameters,pcov)
# #print(f'RC_b = {parameters*10**6} [\mu s]')
# print(f'RC_b = {RC_B*10**6} [\mu s]')
# #plt.plot(np.log(f), A/2.8, 'rx', label='Daten')
# plt.errorbar(np.log(w), A/2.8 , xerr = stds(W_f), yerr = stds(A_f), fmt = 'r.', label='Daten')
# plt.plot(np.log(w), f1(w,*parameters), 'b', label='Fit')

# plt.xlabel(r'$\symup{ln}(\omega)$ [Hz]')
# plt.ylabel(r'$A/U_{0}$')
# plt.legend(loc='best')

# plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
# plt.savefig('build/plot_b.pdf')
# plt.close()

# ## Beginn Plot C

# a_f = unp.uarray(a,0.1) #Fehler von a
# b_f = unp.uarray(b,0.1) #Fehler von b

# phi = a_f/b_f * 2 * np.pi

# def f2(w,c):
#     return np.arctan(-w*c)

# parameters, pcov = curve_fit(f2, w, noms(phi), sigma=None)
# RC_C=unp.uarray(parameters,pcov)
# print(f'RC_c = {-RC_C*10**6} [\mu s]')
# #plt.plot(np.log(f), phi, 'rx', label='Daten')
# plt.errorbar(np.log(w), noms(phi) , xerr = stds(W_f), yerr = stds(phi), fmt = 'r.', label='Daten')
# plt.plot(np.log(w), f2(w,*parameters), 'b', label='Fit')

# plt.xlabel(r'$\symup{log}(\omega)$ [Hz]')
# plt.ylabel(r'$\varphi$ [rad]')
# plt.legend(loc='best')

# plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
# plt.savefig('build/plot_c.pdf')
# plt.close()

# ## Beginn Plot D
# phi2 = np.linspace(0, (np.pi)/2, 50)

# plt.polar(noms(phi), A, 'rx')
# #plt.polar(phi2, np.sin(phi2)/(np.tan(phi2))*2.8, 'b')
# plt.polar(phi2, np.cos(phi2)*2.8, 'b')

# plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
# plt.savefig('build/plot_d.pdf')
# plt.close()
