import matplotlib.pyplot as plt
import numpy as np
import uncertainties as unp
from scipy import stats
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)

#runder Stab einseitig

X_ohne_Fehler, D_0_rund_ohne_Fehler, D_G_rund_ohne_Fehler= np.genfromtxt('content/Daten/einseitig_rund.txt', unpack = True)
X = unp.uarray(X_ohne_Fehler*0.01, 0.002) ###0.001
L = 0.48
Eta = L*X**2-(1/3)*X**3
D_0_rund = unp.uarray(D_0_rund_ohne_Fehler*0.001, 0.00005) ###0.00001
D_G_rund = unp.uarray(D_G_rund_ohne_Fehler*0.001, 0.00005) ###0.00001

D_X = (D_0_rund-D_G_rund)

m_rund_einseitig , b_rund_einseitig , r ,p ,std =stats.linregress(noms(Eta),noms(D_X))

M_rund_einseitig=unp.uarray(m_rund_einseitig,std) #M mit Fehler
B_rund_einseitig=unp.uarray(b_rund_einseitig,std) #B mit Fehler

print(f'M rund einseitig: {M_rund_einseitig}')
print(f'B rund einseitig: {B_rund_einseitig}')
print(f'Eta: {Eta} \n D(x): {D_X}')
X = np.linspace(0, 0.08, 100)

plt.plot(X,m_rund_einseitig*X+b_rund_einseitig ,'b', label = 'Fit')
plt.errorbar(noms(Eta), noms(D_X) , xerr = stds(Eta), yerr = stds(D_X), fmt = 'r.', label='Daten')
plt.xlabel(r'$\eta(x)$ [$\unit{\cubic\meter}$]')
plt.ylabel(r'$D(x)$ [$\unit{\meter}$]')
plt.legend(loc='best')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_rund_einseitig.pdf')
plt.close()

#runder Stab beidseitig

L = 0.54

X1_ohne_Fehler, D1_0_rund_ohne_Fehler, D1_G_rund_ohne_Fehler= np.genfromtxt('content/Daten/beidseitig_rund_1.txt', unpack = True)
X1 = unp.uarray(X1_ohne_Fehler*0.01, 0.001)
Eta_1 = 3*L**2*X1-4*X1**3
D1_0_rund = unp.uarray(D1_0_rund_ohne_Fehler*0.001, 0.00001)
D1_G_rund = unp.uarray(D1_G_rund_ohne_Fehler*0.001, 0.00001)

X2_ohne_Fehler, D2_0_rund_ohne_Fehler, D2_G_rund_ohne_Fehler= np.genfromtxt('content/Daten/beidseitig_rund_2.txt', unpack = True)
X2 = unp.uarray(X2_ohne_Fehler*0.01, 0.001)
Eta_2 = 4*X2**3-12*L*X2**2+9*L**2*X2-L**3
D2_0_rund = unp.uarray(D2_0_rund_ohne_Fehler*0.001, 0.00001)
D2_G_rund = unp.uarray(D2_G_rund_ohne_Fehler*0.001, 0.00001)

D1_X = (D1_0_rund-D1_G_rund)
D2_X = (D2_0_rund-D2_G_rund)

X = np.linspace(0, 0.2, 100)

m1_rund_beidseitig , b1_rund_beidseitig , r ,p ,std =stats.linregress(noms(Eta_1),noms(D1_X))   ####Fit M1

M1_rund_beidseitig=unp.uarray(m1_rund_beidseitig,std) #M1 mit Fehler
B1_rund_beidseitig=unp.uarray(b1_rund_beidseitig,std) #B1 mit Fehler


m2_rund_beidseitig , b2_rund_beidseitig , r ,p ,std =stats.linregress(noms(Eta_2),noms(D2_X))   ####Fit M2

M2_rund_beidseitig=unp.uarray(m2_rund_beidseitig,std) #M2 mit Fehler
B2_rund_beidseitig=unp.uarray(b2_rund_beidseitig,std) #B2 mit Fehler

print(f'M1 rund beidseitig: {M1_rund_beidseitig}')
print(f'M2 rund beidseitig: {M2_rund_beidseitig}')
print(f'B1 rund beidseitig: {B1_rund_beidseitig}')
print(f'B2 rund beidseitig: {B2_rund_beidseitig}')

plt.subplot(1,2,1)
plt.plot(X,m1_rund_beidseitig*X+b1_rund_beidseitig ,'b', label = 'Fit 1')  ####Gerade 1
plt.errorbar(noms(Eta_1), noms(D1_X) , xerr = stds(Eta_1), yerr = stds(D1_X), fmt = 'r.', label='Daten 1')
plt.xlabel(r'$\eta(x)$ [$\unit{\cubic\meter}$]')
plt.ylabel(r'$D(x)$ [$\unit{\meter}$]')
plt.legend(loc='best')
plt.subplot(1,2,2)
plt.plot(X,m2_rund_beidseitig*X+b2_rund_beidseitig ,'m', label = 'Fit 2')  ####Gerade 2
plt.errorbar(noms(Eta_2), noms(D2_X) , xerr = stds(Eta_2), yerr = stds(D2_X), fmt = 'g.', label='Daten 2')
plt.xlabel(r'$\eta(x)$ [$\unit{\cubic\meter}$]')
plt.ylabel(r'$D(x)$ [$\unit{\meter}$]')
plt.legend(loc='best')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_rund_beidseitig.pdf')
plt.close()

#eckiger Stab einseitig

X_ohne_Fehler, D_0_eckig_ohne_Fehler, D_G_eckig_ohne_Fehler= np.genfromtxt('content/Daten/einseitig_eckig.txt', unpack = True)
X = unp.uarray(X_ohne_Fehler*0.01, 0.0001)
L = 0.48
Eta = L*X**2-(1/3)*X**3
D_0_eckig = unp.uarray(D_0_eckig_ohne_Fehler*0.001, 0.00001)
D_G_eckig = unp.uarray(D_G_eckig_ohne_Fehler*0.001, 0.00001)

D_X = (D_0_eckig-D_G_eckig)

m_eckig_einseitig , b_eckig_einseitig , r ,p ,std =stats.linregress(noms(Eta),noms(D_X))

M_eckig_einseitig=unp.uarray(m_eckig_einseitig,std) #M mit Fehler
B_eckig_einseitig=unp.uarray(b_eckig_einseitig,std) #B mit Fehler

print(f'M eckig einseitig: {M_eckig_einseitig}')
print(f'B eckig einseitig: {B_eckig_einseitig}')

X = np.linspace(0, 0.08, 100)

plt.plot(X,m_eckig_einseitig*X+b_eckig_einseitig ,'b', label = 'Fit')
plt.errorbar(noms(Eta), noms(D_X) , xerr = stds(Eta), yerr = stds(D_X), fmt = 'r.', label='Daten')
plt.xlabel(r'$\eta(x)$ [$\unit{\cubic\meter}$]')
plt.ylabel(r'$D(x)$ [$\unit{\meter}$]')
plt.legend(loc='best')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_eckig_einseitig.pdf')
plt.close()

#eckiger Stab beidseitig

L = 0.54

X1_ohne_Fehler, D1_0_eckig_ohne_Fehler, D1_G_eckig_ohne_Fehler= np.genfromtxt('content/Daten/beidseitig_eckig_1.txt', unpack = True)
X1 = unp.uarray(X1_ohne_Fehler*0.01, 0.0001)
Eta_1 = 3*L**2*X1-4*X1**3

D1_0_eckig = unp.uarray(D1_0_eckig_ohne_Fehler*0.001, 0.00001)
D1_G_eckig = unp.uarray(D1_G_eckig_ohne_Fehler*0.001, 0.00001)

X2_ohne_Fehler, D2_0_eckig_ohne_Fehler, D2_G_eckig_ohne_Fehler= np.genfromtxt('content/Daten/beidseitig_eckig_2.txt', unpack = True)
X2 = unp.uarray(X2_ohne_Fehler*0.01, 0.001)
Eta_2 = 4*X2**3-12*L*X2**2+9*L**2*X2-L**3
D2_0_eckig = unp.uarray(D2_0_eckig_ohne_Fehler*0.001, 0.00001)
D2_G_eckig = unp.uarray(D2_G_eckig_ohne_Fehler*0.001, 0.00001)

D1_X = (D1_0_eckig-D1_G_eckig)
D2_X = (D2_0_eckig-D2_G_eckig)

X = np.linspace(0, 0.2, 100)

m1_eckig_beidseitig , b1_eckig_beidseitig , r ,p ,std =stats.linregress(noms(Eta_1),noms(D1_X))   ####Fit M1

M1_eckig_beidseitig=unp.uarray(m1_eckig_beidseitig,std) #M1 mit Fehler
B1_eckig_beidseitig=unp.uarray(b1_eckig_beidseitig,std) #B1 mit Fehler


m2_eckig_beidseitig , b2_eckig_beidseitig , r ,p ,std =stats.linregress(noms(Eta_2),noms(D2_X))   ####Fit M2

M2_eckig_beidseitig=unp.uarray(m2_eckig_beidseitig,std) #M2 mit Fehler
B2_eckig_beidseitig=unp.uarray(b2_eckig_beidseitig,std) #B2 mit Fehler

print(f'M1 eckig beidseitig: {M1_eckig_beidseitig}')
print(f'M2 eckig beidseitig: {M2_eckig_beidseitig}')
print(f'B1 eckig beidseitig: {B1_eckig_beidseitig}')
print(f'B2 eckig beidseitig: {B2_eckig_beidseitig}')

plt.subplot(1,2,1)
plt.plot(X,m1_eckig_beidseitig*X+b1_eckig_beidseitig ,'b', label = 'Fit 1')  ####Gerade 1
plt.errorbar(noms(Eta_1), noms(D1_X) , xerr = stds(Eta_1), yerr = stds(D1_X), fmt = 'r.', label='Daten 1')
plt.xlabel(r'$\eta(x)$ [$\unit{\cubic\meter}$]')
plt.ylabel(r'$D(x)$ [$\unit{\meter}$]')
plt.legend(loc='best')
plt.subplot(1,2,2)
plt.plot(X,m2_eckig_beidseitig*X+b2_eckig_beidseitig ,'m', label = 'Fit 2')  ####Gerade 2
plt.errorbar(noms(Eta_2), noms(D2_X) , xerr = stds(Eta_2), yerr = stds(D2_X), fmt = 'g.', label='Daten 2')
plt.xlabel(r'$\eta(x)$ [$\unit{\cubic\meter}$]')
plt.ylabel(r'$D(x)$ [$\unit{\meter}$]')
plt.legend(loc='best')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_eckig_beidseitig.pdf')
plt.close()

#Elastizitätsmodul loelloeloel#

r = ufloat(0.005, 0.0001)
a = ufloat(0.01, 0.0001)

I_rund = (np.pi/4)*r**4
I_eckig =(1/12)*a**4

F_rund_einseitig = 0.550*9.81
F_rund_beidseitig = 1.750*9.81
F_eckig_einseitig = 0.750*9.81
F_eckig_beidseitig = 1.750*9.81


E_rund_einseitig = F_rund_einseitig/(2*I_rund*M_rund_einseitig)
E1_rund_beidseitig = F_rund_beidseitig/(48*I_rund*M1_rund_beidseitig)
E2_rund_beidseitig = F_rund_beidseitig/(48*I_rund*M2_rund_beidseitig)

E_eckig_einseitig = F_eckig_einseitig/(2*I_eckig*M_eckig_einseitig)
E1_eckig_beidseitig = F_eckig_beidseitig/(48*I_eckig*M1_eckig_beidseitig)
E2_eckig_beidseitig = F_eckig_beidseitig/(48*I_eckig*M2_eckig_beidseitig)

E_rund_gesamt = (E_rund_einseitig + E1_rund_beidseitig + E2_rund_beidseitig)/3
E_eckig_gesamt = (E_eckig_einseitig + E1_eckig_beidseitig + E2_eckig_beidseitig)/3

print('\n')
print('\n')
print(f'Elastizitätsmodul rund einseitig: {E_rund_einseitig}')
print(f'Elastizitätsmodul rund beidseitig 1: {E1_rund_beidseitig}')
print(f'Elastizitätsmodul rund beidseitig 2: {E2_rund_beidseitig}')
print('\n')
print('\n')
print(f'Elastizitätsmodul eckig einseitig: {E_eckig_einseitig}')
print(f'Elastizitätsmodul eckig beidseitig 1: {E1_eckig_beidseitig}')
print(f'Elastizitätsmodul eckig beidseitig 2: {E2_eckig_beidseitig}')
print('\n')
print('\n')
print(f'Elastizitätsmodul rund {E_rund_gesamt}')
print(f'Elastizitätsmodul eckig {E_eckig_gesamt}')
print('\n')
print('\n')