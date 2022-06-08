import matplotlib.pyplot as plt
import numpy as np
import uncertainties as unc
from uncertainties import unumpy as unp 
from uncertainties import ufloat
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as stds
import scipy.optimize as op


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

# %%%%% beta-Strahler %%%%%
# Nullmessung
N_0_beta_unbereinigt = 623
N_0_beta = ufloat(N_0_beta_unbereinigt, np.sqrt(N_0_beta_unbereinigt))/900
d_beta_, d_beta_err, delta_t_beta, N_beta_ = np.genfromtxt('content/data/data_beta.txt', unpack = True)

# Groessen mit Fehlern versehen
d_beta = unp.uarray(d_beta_, d_beta_err)
N_beta = unp.uarray(N_beta_, np.sqrt(N_beta_))

# Normierung von N, korrigieren um Nullmessung ????
N_beta = N_beta / delta_t_beta #- N_0_beta

# N_beta logarithmieren
N_beta_log = unp.log(N_beta)

# cut parameter
cut_beta = 6
d_max = 252

# Plot
plt.errorbar(noms(d_beta), noms(N_beta_log), linestyle = None, fmt='.', xerr = stds(d_beta), yerr = stds(N_beta_log), c='tomato', label='Messwerte')
M_beta, B_beta = regression(noms(d_beta[:cut_beta]),noms(N_beta_log[:cut_beta]),95,270,'dodgerblue','Fit') # Fit
xx = np.linspace(95,485, 1000)
plt.plot(xx, 0*xx+np.log(noms(N_0_beta)), c='slategray', label='Hintergrund')
plt.vlines(d_max, ymin = -0.8, ymax = 3.8, color='darkolivegreen', linestyle = ':', 
label = r'$d_\text{max}=\qty{252}{\micro\metre}$')

plt.xlabel(r'$d \mathbin{/} \unit{\micro\metre}$')
plt.ylabel(r'$\ln(N_{\symup{S}})$')

plt.grid()
plt.legend()
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)

plt.savefig('build/beta.pdf')
plt.close()


# %%%%%% gamma-Strahler %%%%%
# Nullmessung
N_0_gamma_unbereinigt = 1000
N_0_gamma = ufloat(N_0_gamma_unbereinigt, np.sqrt(N_0_gamma_unbereinigt))/900

d_Zn, t_Zn, N_Zn_ = np.genfromtxt('content/data/data_Zn_gamma.txt', unpack = True)
d_Pb, t_Pb, N_Pb_ = np.genfromtxt('content/data/data_Pb_gamma.txt', unpack = True)

# mit Fehlern
N_Zn = unp.uarray(N_Zn_, np.sqrt(N_Zn_))
N_Pb = unp.uarray(N_Pb_, np.sqrt(N_Pb_))

# Normieren, korrigieren um Nullmessung
N_Zn = N_Zn / t_Zn - N_0_gamma
N_Pb = N_Pb / t_Pb - N_0_gamma

# Logarithmieren
N_Zn_log = unp.log(N_Zn)
N_Pb_log = unp.log(N_Pb)

# Plot Zn
plt.errorbar(d_Zn, noms(N_Zn_log), yerr=stds(N_Zn_log), linestyle = None, fmt='.', c='tomato', label='Messwerte')
M_Zn, B_Zn = regression(noms(d_Zn),noms(N_Zn_log),1.5,20.5,'dodgerblue','Fit') # Fit

plt.xlabel(r'$d \mathbin{/} \unit{\milli\metre}$')
plt.ylabel(r'$\ln(N_{\symup{S}})$')

plt.grid()
plt.legend()
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)

plt.savefig('build/Zn.pdf')
plt.close()

# Plot Pb
plt.errorbar(d_Pb, noms(N_Pb_log), yerr=stds(N_Pb_log), linestyle = None, fmt='.', c='tomato', label='Messwerte')
M_Pb, B_Pb = regression(noms(d_Pb),noms(N_Pb_log),-1,41,'dodgerblue','Fit') # Fit

plt.xlabel(r'$d \mathbin{/} \unit{\milli\metre}$')
plt.ylabel(r'$\ln(N_{\symup{S}})$')

plt.grid()
plt.legend()
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)

plt.savefig('build/Pb.pdf')
plt.close()

# PRINT
print('############ Ausgabe V704 ############')
print('------------ Nullmessung -------------')
print(f'N_0_beta: \t',N_0_beta)
print(f'N_0_beta_roh \t',N_0_beta_unbereinigt)
print(f'N_0_gamma \t',N_0_gamma)
print(f'N_0_gamma_roh \t',N_0_gamma_unbereinigt)

print('--------- Plot beta-Strahler ---------')
print(f'M_Zn: \t', M_Zn)
print(f'B_Zn: \t', B_Zn)

print('--------- Plot gamma-Strahler --------')
print(f'M_Pb: \t', M_Pb)
print(f'B_Pb: \t', B_Pb)

