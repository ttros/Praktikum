import matplotlib.pyplot as plt
import numpy as np
import uncertainties as unp
from scipy import stats
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)

#runder Stab einseitig

X_ohne_Fehler, D_0_rund_ohne_Fehler, D_G_rund_ohne_Fehler= np.genfromtxt('content/Daten/einseitig_rund.txt', unpack = True)
X = unp.uarray(X_ohne_Fehler, 0.1)
L = 48
Eta = L*X**2-(1/3)*X**3
D_0_rund = unp.uarray(D_0_rund_ohne_Fehler, 0.01)
D_G_rund = unp.uarray(D_G_rund_ohne_Fehler, 0.01)

D_X = (D_0_rund-D_G_rund)

plt.errorbar(noms(Eta), noms(D_X) , xerr = stds(Eta), yerr = stds(D_X), fmt = 'r.', label='Daten')
plt.xlabel(r'$\eta(x)$ [$\unit{\centi\meter}$]')
plt.ylabel(r'$\Delta D(x)$ [$\unit{\milli\meter}$]')
plt.legend(loc='best')

plt.savefig('build/plot_rund_einseitig.pdf')
plt.close()

#runder Stab beidseitig

X_ohne_Fehler, D_0_rund_ohne_Fehler, D_G_rund_ohne_Fehler= np.genfromtxt('content/Daten/beidseitig_rund.txt', unpack = True)
X = unp.uarray(X_ohne_Fehler, 0.1)
D_0_rund = unp.uarray(D_0_rund_ohne_Fehler, 0.01)
D_G_rund = unp.uarray(D_G_rund_ohne_Fehler, 0.01)

D_X = (D_0_rund-D_G_rund)

plt.errorbar(noms(X), noms(D_X) , xerr = stds(X), yerr = stds(D_X), fmt = 'r.', label='Daten')
plt.xlabel(r'$x$ [$\unit{\centi\meter}$]')
plt.ylabel(r'$\Delta D(x)$ [$\unit{\milli\meter}$]')
plt.legend(loc='best')

plt.savefig('build/plot_rund_beidseitig.pdf')
plt.close()

#eckiger Stab einseitig

X_ohne_Fehler, D_0_eckig_ohne_Fehler, D_G_eckig_ohne_Fehler= np.genfromtxt('content/Daten/einseitig_eckig.txt', unpack = True)
X = unp.uarray(X_ohne_Fehler, 0.1)
L = 48
Eta = L*X**2-(1/3)*X**3
D_0_eckig = unp.uarray(D_0_eckig_ohne_Fehler, 0.01)
D_G_eckig = unp.uarray(D_G_eckig_ohne_Fehler, 0.01)

D_X = (D_0_eckig-D_G_eckig)

plt.errorbar(noms(Eta), noms(D_X) , xerr = stds(Eta), yerr = stds(D_X), fmt = 'r.', label='Daten')
plt.xlabel(r'$\eta(x)$ [$\unit{\centi\meter}$]')
plt.ylabel(r'$\Delta D(x)$ [$\unit{\milli\meter}$]')
plt.legend(loc='best')

plt.savefig('build/plot_eckig_einseitig.pdf')
plt.close()

#eckiger Stab beidseitig

X_ohne_Fehler, D_0_eckig_ohne_Fehler, D_G_eckig_ohne_Fehler= np.genfromtxt('content/Daten/beidseitig_eckig.txt', unpack = True)
X = unp.uarray(X_ohne_Fehler, 0.1)
D_0_eckig = unp.uarray(D_0_eckig_ohne_Fehler, 0.01)
D_G_eckig = unp.uarray(D_G_eckig_ohne_Fehler, 0.01)

D_X = (D_0_eckig-D_G_eckig)

plt.errorbar(noms(X), noms(D_X) , xerr = stds(X), yerr = stds(D_X), fmt = 'r.', label='Daten')
plt.xlabel(r'$x$ [$\unit{\centi\meter}$]')
plt.ylabel(r'$\Delta D(x)$ [$\unit{\milli\meter}$]')
plt.legend(loc='best')

plt.savefig('build/plot_eckig_beidseitig.pdf')
plt.close()