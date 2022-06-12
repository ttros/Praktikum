import matplotlib.pyplot as plt
import numpy as np
import uncertainties as unc
from uncertainties import unumpy as unp 
from uncertainties import ufloat
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as stds
import scipy.optimize as op
import scipy.constants as const


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

# Normierung von N, korrigieren um Nullmessung ????
N_beta_norm = N_beta_ / delta_t_beta #- N_0_beta
print(f'N_beta_norm: {N_beta_norm}')
# Groessen mit Fehlern versehen
d_beta = unp.uarray(d_beta_, d_beta_err)
N_beta = unp.uarray(N_beta_norm, np.sqrt(N_beta_norm))
print(f'N_beta: {N_beta}')

# N_beta logarithmieren
N_beta_log = unp.log(N_beta)
print(f'N_beta_log: {N_beta_log}')
# cut parameter
cut_beta = 6

# Plot
plt.errorbar(noms(d_beta), noms(N_beta_log), linestyle = None, fmt='.', xerr = stds(d_beta), yerr = stds(N_beta_log), c='tomato', capsize=3, label='Messwerte')
M_beta, B_beta = regression(noms(d_beta[:cut_beta]),noms(N_beta_log[:cut_beta]),95,270,'dodgerblue','Fit') # Fit
d_max = (unp.log(N_0_beta) - B_beta) / M_beta # d_max als schnittpunkt berechnen
xx = np.linspace(95,485, 1000)
plt.plot(xx, 0*xx+np.log(noms(N_0_beta)), c='slategray', label='Hintergrund')
plt.vlines(noms(d_max), ymin = -0.8, ymax = 3.8, color='darkolivegreen', linestyle = ':', 
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

# Normieren
N_Zn_norm = N_Zn_ / t_Zn
N_Pb_norm = N_Pb_ / t_Pb

# mit Fehlern, korrigiere um Nullmessung
N_Zn = unp.uarray(N_Zn_norm, np.sqrt(N_Zn_norm))-N_0_gamma
N_Pb = unp.uarray(N_Pb_norm, np.sqrt(N_Pb_norm))-N_0_gamma

# Logarithmieren
N_Zn_log = unp.log(N_Zn)
N_Pb_log = unp.log(N_Pb)

# Plot Zn
plt.errorbar(d_Zn, noms(N_Zn_log), yerr=stds(N_Zn_log), linestyle = None, fmt='.', c='tomato', capsize=3, label='Messwerte')
M_Zn, B_Zn = regression(noms(d_Zn),noms(N_Zn_log),1.5,20.5,'dodgerblue','Fit') # Fit

plt.xlabel(r'$d \mathbin{/} \unit{\milli\metre}$')
plt.ylabel(r'$\ln(N_{\symup{S}})$')

plt.title('Zink')
plt.grid()
plt.legend()
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)

plt.savefig('build/Zn.pdf')
plt.close()

# Plot Pb
plt.errorbar(d_Pb, noms(N_Pb_log), yerr=stds(N_Pb_log), linestyle = None, fmt='.', c='tomato', capsize=3, label='Messwerte')
M_Pb, B_Pb = regression(noms(d_Pb),noms(N_Pb_log),-1,41,'dodgerblue','Fit') # Fit

plt.xlabel(r'$d \mathbin{/} \unit{\milli\metre}$')
plt.ylabel(r'$\ln(N_{\symup{S}})$')

plt.title('Blei')
plt.grid()
plt.legend()
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)

plt.savefig('build/Pb.pdf')
plt.close()

# RECHNUNGEN
r_max = 2.7 * d_max*10**-5
E_max = 1.92 * unp.sqrt((r_max)**2+0.22*r_max)

N_0_Zn = (np.e)**B_Zn
N_0_Pb = (np.e)**B_Pb

# Konstanten

Z_pb = 82
Z_zn = 30
rho_pb = 11300 # kg/m^3
rho_zn =  7140 # kg/m^3
M_pb =  0.2072 # kg/mol
M_zn = 0.06539 # kg/mol

# Theoriewerte
epsilon = 1.295
r_e = 2.82e-15 # m (klassischer Elektronenradius)

sigma_com = 2*np.pi*r_e**2*((1+ epsilon)/epsilon**2 * (2*(1+epsilon)/(1+2*epsilon) - 1/epsilon*np.log(1+2*epsilon)) + 1/(2*epsilon)*np.log(1+2*epsilon) - (1+ 3*epsilon)/(1+2*epsilon)**2)

n_pb = Z_pb*const.N_A*rho_pb/M_pb
n_zn = Z_zn*const.N_A*rho_zn/M_zn

mu_pb_theo = sigma_com*n_pb
mu_zn_theo = sigma_com*n_zn
print("--------------------------------------------------------------------------------")
print("Theoriewerte:")
print("Sigma_com: ", sigma_com, "[m^2]")
print("Blei:    mu_pb = ", mu_pb_theo, "[m^{-1}]")
print("Zink:    mu_zn = ", mu_zn_theo, "[m^{-1}]")
print("--------------------------------------------------------------------------------")

# PRINT
print('############ Ausgabe V704 ############')
print('------------ Nullmessung -------------')
print(f'N_0_beta: \t',N_0_beta)
print(f'N_0_beta_roh \t',N_0_beta_unbereinigt)
print(f'N_0_gamma \t',N_0_gamma)
print(f'N_0_gamma_roh \t',N_0_gamma_unbereinigt)

print('--------- Plot beta-Strahler ---------')
print(f'M_beta: \t', M_beta)
print(f'B_beta: \t', B_beta)
print(f'd_max: \t', d_max)
print(f'r_max: \t', r_max)
print(f'E_max: \t',E_max)
print('--------- Plot gamma-Strahler --------')
print(f'M_Zn: \t', M_Zn)
print(f'B_Zn: \t', B_Zn)
print(f'M_Pb: \t', M_Pb)
print(f'B_Pb: \t', B_Pb)
print(f'N_0_Zn: \t',N_0_Zn)
print(f'N_0_Pb: \t',N_0_Pb)


'''''
############ Ausgabe V704 ############
------------ Nullmessung -------------
N_0_beta: 	 0.692+/-0.028
N_0_beta_roh 	 623
N_0_gamma 	 1.111+/-0.035
N_0_gamma_roh 	 1000
--------- Plot beta-Strahler ---------
M_beta: 	 -0.0234+/-0.0026
B_beta: 	 5.5+/-0.4
--------- Plot gamma-Strahler --------
M_Zn: 	 -0.0426+/-0.0014
B_Zn: 	 4.719+/-0.017
M_Pb: 	 -0.0929+/-0.0019
B_Pb: 	 4.70+/-0.04


############ Ausgabe V704 ############
------------ Nullmessung -------------
N_0_beta: 	 0.692+/-0.028
N_0_beta_roh 	 623
N_0_gamma 	 1.111+/-0.035
N_0_gamma_roh 	 1000
--------- Plot beta-Strahler ---------
M_beta: 	 -0.0234+/-0.0026
B_beta: 	 5.5+/-0.4
d_max: 	 252+/-34
r_max: 	 0.0068+/-0.0009
E_max: 	 0.075+/-0.005
--------- Plot gamma-Strahler --------
M_Zn: 	 -0.0426+/-0.0014
B_Zn: 	 4.719+/-0.017
M_Pb: 	 -0.0929+/-0.0019
B_Pb: 	 4.70+/-0.04
N_0_Zn: 	 112.0+/-2.0
N_0_Pb: 	 110+/-4

'''''