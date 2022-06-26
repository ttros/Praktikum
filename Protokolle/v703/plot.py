import matplotlib.pyplot as plt
import numpy as np
import uncertainties as unc
from uncertainties import unumpy as unp 
from uncertainties import ufloat
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as stds
import scipy.optimize as op
import scipy.constants as const

# lineare Regression
def linfit(x,m,b):
    return m*x+b

def regression(x,y,x1,x2,Farbe,Name):
    params, pcov = op.curve_fit(linfit, x, y)
    std = np.sqrt(np.diag(pcov))
    m = params[0]
    b = params[1]
    M=unp.uarray(m,std[0])
    B=unp.uarray(b,std[1])

    xx = np.linspace(x1,x2,1000)
    plt.plot(xx,m*xx+b,color = Farbe, label=Name)
    return M,B

# %%%%% Charakteristik %%%%%

# Variablen
messdauer = 120   # ggf. normieren, messdauer vermutlich 120s
p_min = 6       # plateau definieren
p_max = -2

# Daten einlesen
U, N_, I_= np.genfromtxt('content/data/data.txt', unpack = True)

N = unp.uarray(N_, np.sqrt(N_)) / messdauer
I = unp.uarray(I_, 0.1)*10**(-6)

# Plot roh
dumme_werte = [10,11,12,13,14,20,21,22,23]
U_mod = np.delete(U,dumme_werte)
N_mod = np.delete(N,dumme_werte)

plt.errorbar(U, noms(N), yerr=stds(N), linestyle = None, fmt='.', c='darkorange', capsize=3, label='Messwerte (Ausreißer)')
plt.errorbar(U_mod, noms(N_mod), yerr=stds(N_mod), linestyle = None, fmt='.', c='forestgreen', capsize=3, label='Messwerte')

plt.xlabel(r'$U \mathbin{/} \unit{\volt}')
plt.ylabel(r'$N \mathbin{/} \unit{\second}$')

plt.grid()
plt.legend()
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)

plt.savefig('build/charakteristik_roh.pdf')
plt.close()

#Plot fein, ohne dummen Werte
plt.errorbar(U_mod[:p_min], noms(N_mod)[:p_min], yerr=stds(N_mod)[:p_min], linestyle = None, fmt='.', c='grey', capsize=3, label='Messwerte')
plt.errorbar(U_mod[p_min:p_max], noms(N_mod)[p_min:p_max], yerr=stds(N_mod)[p_min:p_max], linestyle = None, fmt='.', c='indianred', capsize=3, label='Messwerte Plateaubereich')
plt.errorbar(U_mod[p_max:], noms(N_mod)[p_max:], yerr=stds(N_mod)[p_max:], linestyle = None, fmt='.', c='grey', capsize=3)

M, B = regression(U_mod[p_min:p_max], noms(N_mod)[p_min:p_max], U[p_min]-5, U[p_max-1]+5, 'dodgerblue', 'Fit')

plt.xlabel(r'$U \mathbin{/} \unit{\volt}')
plt.ylabel(r'$N \mathbin{/} \unit{\second}$')

plt.grid()
plt.legend()
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)

plt.savefig('build/charakteristik.pdf')
plt.close()

# %%%%% Totzeit %%%%%
# Daten
N_1_ = 19706
N_2_ = 15897
N_12_ = 34902

# Berechnen
N_1 = ufloat(N_1_, np.sqrt(N_1_))/120
N_2 = ufloat(N_2_, np.sqrt(N_2_))/120
N_12 = ufloat(N_12_, np.sqrt(N_12_))/120

def totzeit(N_1, N_2, N_12):
    return (N_1 + N_2 - N_12)/(2*N_1*N_2)

T = totzeit(N_1, N_2, N_12)

# %%%%% freigesetzte Ladungen %%%%%
dQ = I/(N*const.e)

plt.errorbar(noms(I)*10**6, noms(dQ)*10**-9, yerr=stds(dQ)*10**-9, linestyle = None, fmt='.', c='indianred', capsize=3, label='Ladung')

M, B = regression(noms(I)*10**6, noms(dQ)*10**-9, 0.05, 0.85, 'dodgerblue', 'lineare Ausgleichsgerade')

plt.xlabel(r'$I \mathbin{/} \unit{\micro\ampere}')
plt.ylabel(r'$Q \mathbin{/} \unit{\giga\electronvolt}$')

plt.grid()
plt.legend()
plt.tight_layout()

plt.savefig('build/ladung.pdf')
plt.close()

# PRINT
print('############ Ausgabe V703 ############')
print('----------- Charakteristik -----------')
print(f'lin. Bereich: \t \t {U[p_min]} V bis {U[p_max-1]} V')
print(f'Plateaulaenge: \t \t {U[p_max-1]-U[p_min]} V')
print(f'Steigung: \t \t {M} in ???')
print(f'y_Achsenabschnitt: \t {B} in ???')
print('--------------- Totzeit --------------')
print(f'Z_1: ', '{0:.0f}'.format(ufloat(N_1_, np.sqrt(N_1_))))
print(f'Z_2: ', '{0:.0f}'.format(ufloat(N_2_, np.sqrt(N_2_))))
print(f'Z_12: ', '{0:.0f}'.format(ufloat(N_12_, np.sqrt(N_12_))))
print()
print(f'N_1: {N_1}')
print(f'N_2: {N_2}')
print(f'N_12: {N_12}')
print(f'Totzeit 2 Quellen Meth.: {T} in s')
print('-------------- Ladungen --------------')
print(f'freigesetzte Ladungen in  10^9 eV')
print(f'*******************')
for i in range(len(dQ)):
    print('{0:.1f}'.format(dQ[i]*10**-9))
print(f'*******************')
print(f'*******************')
print(f'N mit Fehler')
for i in range(len(N)):
    print('{0:.0f}'.format(N[i]*messdauer))
print(f'******************')
print(f'Plateausteigung in % / 100 V: {M*100/ufloat(89,np.sqrt(89*120)/120)}')

'''
############ Ausgabe V703 ############
----------- Charakteristik -----------
lin. Bereich: 	 	 390.0 V bis 680.0 V
Plateaulaenge: 	 	 290.0 V
Steigung: 	 	 0.0070+/-0.0019 in ???
y_Achsenabschnitt: 	 85.0+/-1.0 in ???
--------------- Totzeit --------------
Totzeit 2 Quellen Meth.: (5.542+/-0.029)e-05 in s
-------------- Ladungen --------------
freigesetzte Ladungen in eV
*******************
(7.581+/-7.581)e+09
(1.499+/-0.750)e+10
(1.466+/-0.733)e+10
(1.475+/-0.737)e+10
(1.441+/-0.721)e+10
(1.451+/-0.725)e+10
(1.434+/-0.717)e+10
(2.137+/-0.713)e+10
(2.149+/-0.717)e+10
(2.795+/-0.699)e+10
(2.184+/-0.546)e+10
(1.829+/-0.458)e+10
(1.831+/-0.458)e+10
(1.880+/-0.470)e+10
(2.124+/-0.531)e+10
(3.472+/-0.695)e+10
(3.531+/-0.707)e+10
(3.549+/-0.711)e+10
(3.546+/-0.710)e+10
(4.233+/-0.707)e+10
(4.089+/-0.683)e+10
(3.731+/-0.623)e+10
(3.314+/-0.553)e+10
(3.409+/-0.569)e+10
(4.200+/-0.701)e+10
(4.250+/-0.709)e+10
(4.169+/-0.696)e+10
(4.915+/-0.704)e+10
(4.936+/-0.707)e+10
(4.941+/-0.707)e+10
(4.829+/-0.691)e+10
(4.928+/-0.706)e+10
(4.872+/-0.698)e+10
(5.508+/-0.691)e+10
(5.561+/-0.697)e+10
(5.586+/-0.700)e+10
(5.445+/-0.683)e+10
(5.369+/-0.673)e+10
*******************
*******************
N mit Fehler
9880+/-99
9994+/-100
10221+/-101
10159+/-101
10393+/-102
10326+/-102
10443+/-102
10516+/-103
10457+/-102
10719+/-104
13720+/-117
16377+/-128
16365+/-128
15932+/-126
14108+/-119
10787+/-104
10606+/-103
10553+/-103
10562+/-103
10617+/-103
10990+/-105
12046+/-110
13559+/-116
13181+/-115
10699+/-103
10575+/-103
10779+/-104
10668+/-103
10622+/-103
10611+/-103
10858+/-104
10640+/-103
10761+/-104
10878+/-104
10774+/-104
10727+/-104
11005+/-105
11160+/-106
******************
'''
