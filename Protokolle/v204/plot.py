import matplotlib.pyplot as plt
import numpy as np
import uncertainties as unp
from scipy import stats
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
######
#Markersize
mksz = 2.0
#####

######
# Celsius in Kelvin
def cToK(temp):
    return temp + 273.15
######
######
# statisch
# 1,4,5,8 in einem
t_roh, T_1_c, T_2_c, T_3_c, T_4_c, T_5_c, T_6_c, T_7_c, T_8_c = np.genfromtxt('content/data/data_statisch.txt', unpack = True)
t = t_roh * 5
T_1 = cToK(T_1_c)
T_2 = cToK(T_2_c)
T_3 = cToK(T_3_c)
T_4 = cToK(T_4_c)
T_5 = cToK(T_5_c)
T_6 = cToK(T_6_c)
T_7 = cToK(T_7_c)
T_8 = cToK(T_8_c)

plt.plot(t, T_1, '.', markersize = mksz, label = r'$T_{1}$, Messing, breit')
plt.plot(t, T_4, '.', markersize = mksz, label = r'$T_{4}$, Messing, schmal')
plt.plot(t, T_5, '.', markersize = mksz, label = r'$T_{5}$, Aluminium')
plt.plot(t, T_8, '.', markersize = mksz, label = r'$T_{8}$, Edelstahl')

plt.xlabel(r'$t / \unit{\second}$')
plt.ylabel(r'$T / \unit{\kelvin}$')
plt.grid()
plt.legend(loc='best')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_statisch.pdf')
plt.close()

# Berechnung des Waermestroms
#Waermeleitfaehigkeit kappa
k_me = 88
k_al = 234
k_es = 25

#Querschnitte
a_breit = 1.2 * 0.4 / 10000
a_schmal = 0.7 * 0.4 / 10000

#Abstand x
x = 0.03

def diffT_meb(t):
    return T_2[t]-T_1[t]

T_meb_100 = diffT_meb(20)
T_meb_200 = diffT_meb(40)
T_meb_300 = diffT_meb(60)
T_meb_400 = diffT_meb(80)
T_meb_700 = diffT_meb(140)

Q_meb_100 = - k_me * a_breit * T_meb_100 / x
Q_meb_200 = - k_me * a_breit * T_meb_200 / x
Q_meb_300 = - k_me * a_breit * T_meb_300 / x
Q_meb_400 = - k_me * a_breit * T_meb_400 / x
Q_meb_700 = - k_me * a_breit * T_meb_700 / x

print('\n')
print(f'T_meb_100: {T_meb_100}')
print(f'T_meb_200: {T_meb_200}')
print(f'T_meb_300: {T_meb_300}')
print(f'T_meb_400: {T_meb_400}')
print(f'T_meb_700: {T_meb_700}')
print(f'Q_meb: {Q_meb_100, Q_meb_200, Q_meb_300, Q_meb_400, Q_meb_700}')
print('\n')

def diffT_mes(t):
    return T_3[t]-T_4[t]
T_mes_100 = diffT_mes(20)
T_mes_200 = diffT_mes(40)
T_mes_300 = diffT_mes(60)
T_mes_400 = diffT_mes(80)
T_mes_700 = diffT_mes(140)
Q_mes_100 = - k_me * a_schmal * T_mes_100 / x
Q_mes_200 = - k_me * a_schmal * T_mes_200 / x
Q_mes_300 = - k_me * a_schmal * T_mes_300 / x
Q_mes_400 = - k_me * a_schmal * T_mes_400 / x
Q_mes_700 = - k_me * a_schmal * T_mes_700 / x
print(f'T_mes_100: {T_mes_100}')
print(f'T_mes_200: {T_mes_200}')
print(f'T_mes_300: {T_mes_300}')
print(f'T_mes_400: {T_mes_400}')
print(f'T_mes_700: {T_mes_700}')
print(f'Q_mes: {Q_mes_100, Q_mes_200, Q_mes_300, Q_mes_400, Q_mes_700}')
print('\n')

def diffT_al(t):
    return T_6[t]-T_5[t]
T_al_100 = diffT_al(20)
T_al_200 = diffT_al(40)
T_al_300 = diffT_al(60)
T_al_400 = diffT_al(80)
T_al_700 = diffT_al(140)
Q_al_100 = - k_al * a_breit * T_al_100 / x
Q_al_200 = - k_al * a_breit * T_al_200 / x
Q_al_300 = - k_al * a_breit * T_al_300 / x
Q_al_400 = - k_al * a_breit * T_al_400 / x
Q_al_700 = - k_al * a_breit * T_al_700 / x
print(f'T_al_100: {T_al_100}')
print(f'T_al_200: {T_al_200}')
print(f'T_al_300: {T_al_300}')
print(f'T_al_400: {T_al_400}')
print(f'T_al_700: {T_al_700}')
print(f'Q_al: {Q_al_100, Q_al_200, Q_al_300, Q_al_400, Q_al_700}')
print('\n')

def diffT_es(t):
    return T_7[t]-T_8[t]
T_es_100 = diffT_es(20)
T_es_200 = diffT_es(40)
T_es_300 = diffT_es(60)
T_es_400 = diffT_es(80)
T_es_700 = diffT_es(140)
Q_es_100 = - k_es * a_breit * T_es_100 / x
Q_es_200 = - k_es * a_breit * T_es_200 / x
Q_es_300 = - k_es * a_breit * T_es_300 / x
Q_es_400 = - k_es * a_breit * T_es_400 / x
Q_es_700 = - k_es * a_breit * T_es_700 / x
print(f'T_es_100: {T_es_100}')
print(f'T_es_200: {T_es_200}')
print(f'T_es_300: {T_es_300}')
print(f'T_es_400: {T_es_400}')
print(f'T_es_700: {T_es_700}')
print(f'Q_es: {Q_es_100, Q_es_200, Q_es_300, Q_es_400, Q_es_700}')
print('\n')

#T_7-T_8
T_7_8 = T_7 - T_8
T_2_1 = T_2 - T_1

plt.plot(t, T_7_8, '.', markersize = mksz, label = r'$\symup{\Delta}T_{\symup{Edelstahl}} = T_{7} - T_{8}$')
plt.plot(t, T_2_1, '.', markersize = mksz, label = r'$\symup{\Delta}T_{\symup{Messing}} = T_{2} - T_{1}$')

plt.xlabel(r'$t / \unit{\second}$')
plt.ylabel(r'$\symup{\Delta}T / \unit{\kelvin}$')
plt.grid()
plt.legend(loc='best')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_T_7-8.pdf')
plt.close()

######
# dynamisch, kurz
t_roh, T_1_c, T_2_c, T_3_c, T_4_c, T_5_c, T_6_c, T_7_c, T_8_c = np.genfromtxt('content/data/data_dynamisch_kurz.txt', unpack = True)
t = t_roh * 2
T_1 = cToK(T_1_c)
T_2 = cToK(T_2_c)
T_3 = cToK(T_3_c)
T_4 = cToK(T_4_c)
T_5 = cToK(T_5_c)
T_6 = cToK(T_6_c)
T_7 = cToK(T_7_c)
T_8 = cToK(T_8_c)

plt.plot(t, T_1, '.',markersize = mksz, label = r'$T_{1}$, Messing, breit, fern')
plt.plot(t, T_2, '.', markersize = mksz, label = r'$T_{2}$, Messing, breit, nah')

plt.xlabel(r'$t / \unit{\second}$')
plt.ylabel(r'$T / \unit{\kelvin}$')
plt.grid()
plt.legend(loc='best')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_dynamisch_kurz.pdf')
plt.close()
######

######
# dynamisch, lang
t_roh, T_1_c, T_2_c, T_3_c, T_4_c, T_5_c, T_6_c, T_7_c, T_8_c = np.genfromtxt('content/data/data_dynamisch_lang.txt', unpack = True)
t = t_roh * 2
T_1 = cToK(T_1_c)
T_2 = cToK(T_2_c)
T_3 = cToK(T_3_c)
T_4 = cToK(T_4_c)
T_5 = cToK(T_5_c)
T_6 = cToK(T_6_c)
T_7 = cToK(T_7_c)
T_8 = cToK(T_8_c)

plt.plot(t, T_7, '.',markersize = mksz, label = r'$T_{7}$, Edelstahl, nah')
plt.plot(t, T_8, '.', markersize = 1.0, label = r'$T_{8}$, Edelstahl, fern')

plt.xlabel(r'$t / \unit{\second}$')
plt.ylabel(r'$T / \unit{\kelvin}$')
plt.grid()
plt.legend(loc='best')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_dynamisch_lang.pdf')
plt.close()
######