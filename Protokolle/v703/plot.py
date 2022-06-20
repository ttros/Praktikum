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
p_min = 3       # plateau definieren
p_max = 31

# Daten einlesen
U, N_, I_= np.genfromtxt('content/data/data_muster.txt', unpack = True)

N = unp.uarray(N_, np.sqrt(N_)) / messdauer
I = unp.uarray(I_, 0.1)*10**(-6)

plt.errorbar(U[:p_min], noms(N)[:p_min], yerr=stds(N)[:p_min], linestyle = None, fmt='.', c='grey', capsize=3, label='Messwerte')
plt.errorbar(U[p_min:p_max], noms(N)[p_min:p_max], yerr=stds(N)[p_min:p_max], linestyle = None, fmt='.', c='indianred', capsize=3, label='Messwerte Plateaubereich')
plt.errorbar(U[p_max:], noms(N)[p_max:], yerr=stds(N)[p_max:], linestyle = None, fmt='.', c='grey', capsize=3)

M, B = regression(U[p_min:p_max], noms(N)[p_min:p_max], U[p_min]-5, U[p_max-1]+5, 'dodgerblue', 'Fit')

plt.xlabel(r'$U \mathbin{/} \unit{\volt}')
plt.ylabel(r'$N \mathbin{/} \unit{\second}$')

plt.grid()
plt.legend()
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)

plt.savefig('build/charakteristik.pdf')
plt.close()

# %%%%% Totzeit %%%%%
# Daten
N_1_ = 111
N_2_ = 111
N_12_ = 222

# Berechnen
N_1 = ufloat(N_1_, np.sqrt(N_1_))
N_2 = ufloat(N_2_, np.sqrt(N_2_))
N_12 = ufloat(N_12_, np.sqrt(N_12_))

def totzeit(N_1, N_2, N_12):
    return (N_1 + N_2 - N_12)/(2*N_1*N_2)

T = totzeit(N_1, N_2, N_12)

# %%%%% freigesetzte Ladungen %%%%%
dQ = I/(N*const.e)

plt.errorbar(U, noms(dQ)*10**9, yerr=stds(dQ)*10**9, linestyle = None, fmt='.', c='indianred', capsize=3, label='Ladung')

plt.xlabel(r'$U \mathbin{/} \unit{\volt}')
plt.ylabel(r'$\symup{d}Q \mathbin{/} \unit{\nano\electronvolt}$')

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
print(f'Totzeit 2 Quellen Meth.: {T} in s')
print('-------------- Ladungen --------------')
print(f'freigesetzte Ladungen in eV')
print(f'*******************')
for i in range(len(dQ)):
    print('{0:.3e}'.format(dQ[i]))
print(f'*******************')
