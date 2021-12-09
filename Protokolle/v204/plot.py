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