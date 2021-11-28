import matplotlib.pyplot as plt
import numpy as np
import uncertainties as unp
from scipy import stats
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)

# Apperaturkonstante, Viskosität, Reinoldszahl

# Fallzeiten
# kl. Kugel
O_k_, U_k_ = np.genfromtxt('content/data_kl_Kugel.txt', unpack = True)
O_k = unp.uarray(O_k_, 0.5)
O_k_mittel = sum(O_k)/len(O_k)

U_k = unp.uarray(U_k_, 0.5)
U_k_mittel = sum(U_k)/len(U_k)

print(f'O_k_mean: {O_k_mittel}')
print(f'U_k_mean: {U_k_mittel}')

# gr. Kugel
O_g_, U_g_ = np.genfromtxt('content/data_gr_Kugel.txt', unpack = True)
O_g = unp.uarray(O_g_, 0.5)
O_g_mittel = sum(O_g)/len(O_g)

U_g = unp.uarray(U_g_, 0.5)
U_g_mittel = sum(U_g)/len(U_g)

print(f'O_g_mean: {O_g_mittel}')
print(f'U_g_mean: {U_g_mittel}')

# Dichte der Kugeln
d_k = ufloat(1.565, 0.01) #cm
m_k = 4.4531 #g
dichte_k = m_k / (4/3 * np.pi * ((d_k/2))**3)

d_g = ufloat(1.585, 0.01) #cm
m_g = 4.9528 #g
dichte_g = m_g / (4/3 * np.pi * ((d_g/2))**3)

print(f'Dichte kl. Kugel {dichte_k} in g/cm^3')
print(f'Dichte gr. Kugel {dichte_g} in g/cm^3')

# Viskosität
K = 0.07640
visk_O_k = K*(dichte_k - 0.99823)*2*O_k_mittel
visk_U_k = K*(dichte_k - 0.99823)*2*U_k_mittel

print(f'Viskositaet kl. Kugel oben: {visk_O_k}')
print(f'Viskositaet kl. Kugel unten: {visk_U_k}')

# Apperaturkonstante gr. Kugel
K_gr_O = visk_O_k/((dichte_g - 0.99823)*2*O_g_mittel)
K_gr_U = visk_O_k/((dichte_g - 0.99823)*2*U_g_mittel)

print(f'Apperaturkonstante gr. Kugel oben: {K_gr_O}')
print(f'Apperaturkonstante gr. Kugel unten: {K_gr_U}')

# Reinoldszahl
dichte_w = 0.99823
Re_kl_O = 100 * dichte_w * (5/O_k_mittel) * d_g / visk_O_k
Re_kl_U = 100 * dichte_w * (5/U_k_mittel) * d_g / visk_U_k

Re_gr_O = 100 * dichte_w * (5/O_g_mittel) * d_g / visk_O_k
Re_gr_U = 100 * dichte_w * (5/U_g_mittel) * d_g / visk_U_k

print(f'Reinoldszahl kl. Kugel oben: {Re_kl_O}')
print(f'Reinoldszahl kl. Kugel unten: {Re_kl_U}')
print(f'Reinoldszahl gr. Kugel oben: {Re_gr_O}')
print(f'Reinoldszahl gr. Kugel unten: {Re_gr_U}')



# Plots
# dynamische Viskosität

T_, t_O_1_, t_O_2_, t_U_1_, t_U_2_ = np.genfromtxt('content/data_T.txt', unpack = True)
# Messunsicherheiten
T = unp.uarray(T_,1)

t_O_=(t_O_1_+t_O_2_)/2
t_O = unp.uarray(t_O_,0.5)
t_U_=(t_U_1_+t_U_2_)/2
t_U = unp.uarray(t_U_,0.5)


# Plot Oben
m , b , r ,p ,std =stats.linregress(1/noms(T),np.log(noms(t_O)))
M=unp.uarray(m,std)
B=unp.uarray(b,std)
plt.plot(1/noms(T), m*1/noms(T)+b, 'b', label = 'Fit')
plt.errorbar(1/noms(T), np.log(noms(t_O)), xerr = stds(1/T), yerr = stds(unp.log(t_O)), fmt = 'r.', label='Daten')
plt.xlabel(r'$\frac{1}{T}$ [$\unit{\per\celsius}$]')
plt.ylabel(r'$\symup{ln}(t)$ [$\unit{\second}$]')
plt.legend(loc='best')

print(f'Geradegleichug Oben: {M}*1/T + {B}')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_oben.pdf')
plt.close()

# Plot Unten
m , b , r ,p ,std =stats.linregress(1/noms(T),np.log(noms(t_U)))
M=unp.uarray(m,std)
B=unp.uarray(b,std)
plt.plot(1/noms(T), m*1/noms(T)+b, 'b', label = 'Fit')
plt.errorbar(1/noms(T), np.log(noms(t_U)), xerr = stds(1/T), yerr = stds(unp.log(t_U)), fmt = 'r.', label='Daten')
plt.xlabel(r'$\frac{1}{T}$ [$\unit{\per\celsius}$]')
plt.ylabel(r'$\symup{ln}(t)$ [$\unit{\second}$]')
plt.legend(loc='best')

print(f'Geradegleichug Unten: {M}*1/T + {B}')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_unten.pdf')
plt.close()
