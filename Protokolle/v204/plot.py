import matplotlib.pyplot as plt
import numpy as np
import uncertainties as unp
from scipy import stats
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)

######
# statisch

t_roh, T_1_roh, T_2_roh, T_3_roh, T_4_roh, T_5_roh, T_6_roh, T_7_roh, T_8_roh = np.genfromtxt('content/data/data_statisch.txt', unpack = True)

plt.plot(t_roh, T_1_roh)


#plt.plot(X,m_rund_einseitig*X+b_rund_einseitig ,'b', label = 'Fit')
#plt.errorbar(noms(Eta), noms(D_X) , xerr = stds(Eta), yerr = stds(D_X), fmt = 'r.', label='Daten')
#plt.xlabel(r'$\eta(x)$ [$\unit{\cubic\meter}$]')
#plt.ylabel(r'$D(x)$ [$\unit{\meter}$]')
#plt.legend(loc='best')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_statisch.pdf')
plt.close()

######
# dynamisch, kurz
t_roh, T_1_roh, T_2_roh, T_3_roh, T_4_roh, T_5_roh, T_6_roh, T_7_roh, T_8_roh = np.genfromtxt('content/data/data_dynamisch_kurz.txt', unpack = True)

plt.plot(t_roh, T_1_roh)
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_dynamisch_kurz.pdf')
plt.close()
######

######
t_roh, T_1_roh, T_2_roh, T_3_roh, T_4_roh, T_5_roh, T_6_roh, T_7_roh, T_8_roh = np.genfromtxt('content/data/data_dynamisch_lang.txt', unpack = True)

plt.plot(t_roh, T_1_roh)
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_dynamisch_lang.pdf')
plt.close()
######