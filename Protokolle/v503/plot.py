import matplotlib.pyplot as plt
import numpy as np

import uncertainties as unp
from scipy import stats
from scipy.stats import sem
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
import scipy.constants as const
import scipy.optimize as op

# %%% Definition Konstanten %%%
d = ufloat(7.6250,0.0051)*10**-3            # in m
dichte_oel = 886                            # in kg/m^3
dichte_L = 1.204                            # in kg/m^3 bei 20°C
B = 6.17*10**(-5)*133.322                   # in Pa*m

# %%% Temperaturen in °C und unkorrigierten Viskositäten %%%
T_1 = 27.0
T_2 = 27.5
T_3 = 27.5
T_4 = 28.0
T_5 = 28.0
viskosität = {
    'eta_L_1': 1.855e-5,
    'eta_L_2': 1.860e-5, 
    'eta_L_3': 1.860e-5, 
    'eta_L_4': 1.862e-5, 
    'eta_L_5': 1.862e-5 
}
spannungen = {
    'U_1': 150.1,
    'U_2': 175.4,
    'U_3': 200.9,
    'U_4': 225.0,
    'U_5': 250.0
}

# %%% Definition Formeln %%%

def ladung(v_auf, v_ab, U, eta_L):
    return 9/2*np.pi*unp.sqrt((eta_L**3*abs(v_ab-v_auf))/(const.g*(dichte_oel-dichte_L))) * d*(v_ab + v_auf)/U

def radius(v_auf, v_ab, eta_L):
    return unp.sqrt((9*eta_L*abs(v_ab-v_auf))/(2*(const.g)*(dichte_oel-dichte_L)))

def ladung_korrigiert(v_auf, v_ab, U, eta_L):
    return ladung(v_auf, v_ab, U, eta_L)*(1+B*1/(100000*radius(v_auf, v_ab, eta_L)))

# %%% Bestimmung Ladung der Tröpfchen %%%
def auswertung(messung_nummer):
    t_auf_1, t_ab_1, t_auf_2, t_ab_2, t_auf_3, t_ab_3, t_auf_4, t_ab_4, t_auf_5, t_ab_5 = np.genfromtxt(f'content/data/Messung_{messung_nummer}.txt', unpack = True)

    times = {
        't_auf_1': t_auf_1,
        't_ab_1' : t_ab_1 ,
        't_auf_2': t_auf_2,
        't_ab_2' : t_ab_2 ,
        't_auf_3': t_auf_3,
        't_ab_3' : t_ab_3 ,
        't_auf_4': t_auf_4,
        't_ab_4' : t_ab_4 ,
        't_auf_5': t_auf_5,
        't_ab_5' : t_ab_5 
    }
    velocities = {
        'v_auf_1': 0.5*10**-3/t_auf_1,
        'v_ab_1' : 0.5*10**-3/t_ab_1 ,
        'v_auf_2': 0.5*10**-3/t_auf_2,
        'v_ab_2' : 0.5*10**-3/t_ab_2 ,
        'v_auf_3': 0.5*10**-3/t_auf_3,
        'v_ab_3' : 0.5*10**-3/t_ab_3 ,
        'v_auf_4': 0.5*10**-3/t_auf_4,
        'v_ab_4' : 0.5*10**-3/t_ab_4 ,
        'v_auf_5': 0.5*10**-3/t_auf_5,
        'v_ab_5' : 0.5*10**-3/t_ab_5 
    }

    velocities_mean = {
        'v_m_auf_1': ufloat(np.mean(velocities['v_auf_1']), sem(velocities['v_auf_1'])),
        'v_m_ab_1' : ufloat(np.mean(velocities['v_ab_1' ]), sem(velocities['v_ab_1' ])),
        'v_m_auf_2': ufloat(np.mean(velocities['v_auf_2']), sem(velocities['v_auf_2'])),
        'v_m_ab_2' : ufloat(np.mean(velocities['v_ab_2' ]), sem(velocities['v_ab_2' ])),
        'v_m_auf_3': ufloat(np.mean(velocities['v_auf_3']), sem(velocities['v_auf_3'])),
        'v_m_ab_3' : ufloat(np.mean(velocities['v_ab_3' ]), sem(velocities['v_ab_3' ])),
        'v_m_auf_4': ufloat(np.mean(velocities['v_auf_4']), sem(velocities['v_auf_4'])),
        'v_m_ab_4' : ufloat(np.mean(velocities['v_ab_4' ]), sem(velocities['v_ab_4' ])),
        'v_m_auf_5': ufloat(np.mean(velocities['v_auf_5']), sem(velocities['v_auf_5'])),
        'v_m_ab_5' : ufloat(np.mean(velocities['v_ab_5' ]), sem(velocities['v_ab_5' ]))
    }

    radii = {
        'r_1': radius(velocities_mean['v_m_auf_1'], velocities_mean['v_m_ab_1'], viskosität[f'eta_L_{messung_nummer}']),
        'r_2': radius(velocities_mean['v_m_auf_2'], velocities_mean['v_m_ab_2'], viskosität[f'eta_L_{messung_nummer}']),
        'r_3': radius(velocities_mean['v_m_auf_3'], velocities_mean['v_m_ab_3'], viskosität[f'eta_L_{messung_nummer}']),
        'r_4': radius(velocities_mean['v_m_auf_4'], velocities_mean['v_m_ab_4'], viskosität[f'eta_L_{messung_nummer}']),
        'r_5': radius(velocities_mean['v_m_auf_5'], velocities_mean['v_m_ab_5'], viskosität[f'eta_L_{messung_nummer}'])
    }
    charges = {
        'q_1': ladung(velocities_mean['v_m_auf_1'], velocities_mean['v_m_ab_1'], spannungen[f'U_{messung_nummer}'], viskosität[f'eta_L_{messung_nummer}']),
        'q_2': ladung(velocities_mean['v_m_auf_2'], velocities_mean['v_m_ab_2'], spannungen[f'U_{messung_nummer}'], viskosität[f'eta_L_{messung_nummer}']),
        'q_3': ladung(velocities_mean['v_m_auf_3'], velocities_mean['v_m_ab_3'], spannungen[f'U_{messung_nummer}'], viskosität[f'eta_L_{messung_nummer}']),
        'q_4': ladung(velocities_mean['v_m_auf_4'], velocities_mean['v_m_ab_4'], spannungen[f'U_{messung_nummer}'], viskosität[f'eta_L_{messung_nummer}']),
        'q_5': ladung(velocities_mean['v_m_auf_5'], velocities_mean['v_m_ab_5'], spannungen[f'U_{messung_nummer}'], viskosität[f'eta_L_{messung_nummer}'])
    }
    charges_real = {
        'q_real_1': ladung_korrigiert(velocities_mean['v_m_auf_1'], velocities_mean['v_m_ab_1'], spannungen[f'U_{messung_nummer}'], viskosität[f'eta_L_{messung_nummer}']),
        'q_real_2': ladung_korrigiert(velocities_mean['v_m_auf_2'], velocities_mean['v_m_ab_2'], spannungen[f'U_{messung_nummer}'], viskosität[f'eta_L_{messung_nummer}']),
        'q_real_3': ladung_korrigiert(velocities_mean['v_m_auf_3'], velocities_mean['v_m_ab_3'], spannungen[f'U_{messung_nummer}'], viskosität[f'eta_L_{messung_nummer}']),
        'q_real_4': ladung_korrigiert(velocities_mean['v_m_auf_4'], velocities_mean['v_m_ab_4'], spannungen[f'U_{messung_nummer}'], viskosität[f'eta_L_{messung_nummer}']),
        'q_real_5': ladung_korrigiert(velocities_mean['v_m_auf_5'], velocities_mean['v_m_ab_5'], spannungen[f'U_{messung_nummer}'], viskosität[f'eta_L_{messung_nummer}'])
    }


    # %%% Ausgabe %%%
    print(f'-------------------------------------------------')
    print(f'----------------- Messung Nr. {messung_nummer} -----------------')
    print(f'------------- Geschwindigkeiten -----------------')
    for x in velocities_mean:
        print({x},':\t', '{0:.3e}'.format(velocities_mean[f'{x}']))
    print(f'-------------------- Radien ---------------------')
    for x in radii:
        print({x},':\t', '{0:.3e}'.format(radii[f'{x}']))
    print(f'------------------- Ladungen --------------------')
    for x in charges:
        print({x},':\t', '{0:.3e}'.format(charges[f'{x}']))
    print(f'-------------- Ladungen korrigiert --------------')
    for x in charges_real:
        print({x},':\t', '{0:.3e}'.format(charges_real[f'{x}']))
    print(f'#################################################')

    return charges, charges_real
# %%% Ende Auswertung %%%

ladungen_1_dic, ladungen_1_korr_dic = auswertung(1)
ladungen_2_dic, ladungen_2_korr_dic = auswertung(2)
ladungen_3_dic, ladungen_3_korr_dic = auswertung(3)
ladungen_4_dic, ladungen_4_korr_dic = auswertung(4)
ladungen_5_dic, ladungen_5_korr_dic = auswertung(5)

ladungen_1 = unp.uarray([0,0,0,0,0],[0,0,0,0,0])
ladungen_2 = unp.uarray([0,0,0,0,0],[0,0,0,0,0])
ladungen_3 = unp.uarray([0,0,0,0,0],[0,0,0,0,0])
ladungen_4 = unp.uarray([0,0,0,0,0],[0,0,0,0,0])
ladungen_5 = unp.uarray([0,0,0,0,0],[0,0,0,0,0])

ladungen_1_korr = unp.uarray([0,0,0,0,0],[0,0,0,0,0])
ladungen_2_korr = unp.uarray([0,0,0,0,0],[0,0,0,0,0])
ladungen_3_korr = unp.uarray([0,0,0,0,0],[0,0,0,0,0])
ladungen_4_korr = unp.uarray([0,0,0,0,0],[0,0,0,0,0])
ladungen_5_korr = unp.uarray([0,0,0,0,0],[0,0,0,0,0])


for i in range(0,5):
    ladungen_1[i] = ladungen_1_dic[f'q_{i+1}']
    ladungen_2[i] = ladungen_2_dic[f'q_{i+1}']
    ladungen_3[i] = ladungen_3_dic[f'q_{i+1}']
    ladungen_4[i] = ladungen_4_dic[f'q_{i+1}']
    ladungen_5[i] = ladungen_5_dic[f'q_{i+1}']

    ladungen_1_korr[i] = ladungen_1_korr_dic[f'q_real_{i+1}']
    ladungen_2_korr[i] = ladungen_2_korr_dic[f'q_real_{i+1}']
    ladungen_3_korr[i] = ladungen_3_korr_dic[f'q_real_{i+1}']
    ladungen_4_korr[i] = ladungen_4_korr_dic[f'q_real_{i+1}']
    ladungen_5_korr[i] = ladungen_5_korr_dic[f'q_real_{i+1}']

# Funktion zur Bestimmung von e_0
def elementar(arr):
    e_0 = 999
    for i in range(len(arr)):
        for j in (np.array(range(len(ladungen_1)-1)) + 1):
            a = np.abs(arr[i] - arr[j])
            if a > 1e-19 and a < e_0:
                e_0 = a
    return e_0 

e_0_1 = elementar(ladungen_1)
e_0_2 = elementar(ladungen_2)
e_0_3 = elementar(ladungen_3)
e_0_4 = elementar(ladungen_4)
e_0_5 = elementar(ladungen_5)
# e_0_mean = sum([e_0_1, e_0_2, e_0_3, e_0_4, e_0_5])/5
e_0_mean = (e_0_1 + e_0_2 + e_0_3 + e_0_5)/4

e_0_1_korr = elementar(ladungen_1_korr)
e_0_2_korr = elementar(ladungen_2_korr)
e_0_3_korr = elementar(ladungen_3_korr)
e_0_4_korr = elementar(ladungen_4_korr)
e_0_5_korr = elementar(ladungen_5_korr)
# e_0_korr_mean = sum([e_0_1_korr, e_0_2_korr, e_0_3_korr, e_0_4_korr, e_0_5_korr])/5
e_0_korr_mean = (e_0_1_korr  + e_0_3_korr  + e_0_5_korr)/3

# ladungen_2_korr = np.delete(ladungen_2,1)   # löschen von mülldaten

print(f'-------------------------------------------------')
print(f'--------- gemittelte Elementarladungen ----------')
print(f'e_0_1: \t\t', '{0:.3e}'.format(e_0_1))
print(f'e_0_2: \t\t', '{0:.3e}'.format(e_0_2))
print(f'e_0_3: \t\t', '{0:.3e}'.format(e_0_3))
print(f'e_0_4: \t\t', '{0:.3e}'.format(e_0_4))
print(f'e_0_5: \t\t', '{0:.3e}'.format(e_0_5))
print(f'e_0_mean: \t', '{0:.3e}'.format(e_0_mean))
print(f'-------------------------------------------------')
print(f'----- gemittelte Elementarladung, korrigiert ----')
print(f'e_0_1_korr: \t', '{0:.3e}'.format(e_0_1_korr))
print(f'e_0_2_korr: \t', '{0:.3e}'.format(e_0_2_korr))
print(f'e_0_3_korr: \t', '{0:.3e}'.format(e_0_3_korr))
print(f'e_0_4_korr: \t', '{0:.3e}'.format(e_0_4_korr))
print(f'e_0_5_korr: \t', '{0:.3e}'.format(e_0_5_korr))
print(f'e_0_korr_mean: \t', '{0:.3e}'.format(e_0_korr_mean))
print(f'-------------------------------------------------')


F = ufloat(96485.3399, 0.0024)
N_A = F/e_0_korr_mean
dN_A = np.abs(noms(N_A)- const.N_A)/(const.N_A)
print("Avogadrokonstante: ", '{0:.4e}'.format(N_A))
print("Abweichung: ", '{0:.4f}'.format(dN_A))

# Plot zur Ladungsverteilung
e = const.e
x = np.array([0.9, 0.95, 1, 1.05, 1.1])
# x_korr = np.array([0.9, 0.95, 1, 1.05])
plt.errorbar(x,   noms(ladungen_1), yerr = stds(ladungen_1), elinewidth = 0.7, linewidth = 0, marker = ".", markersize = 7, capsize=3)
# plt.errorbar(1+x_korr, noms(ladungen_2_korr), yerr = stds(ladungen_2_korr), elinewidth = 0.7, linewidth = 0, marker = ".", markersize = 7, capsize=3)
plt.errorbar(1+x, noms(ladungen_2), yerr = stds(ladungen_2), elinewidth = 0.7, linewidth = 0, marker = ".", markersize = 7, capsize=3)
plt.errorbar(2+x, noms(ladungen_3), yerr = stds(ladungen_3), elinewidth = 0.7, linewidth = 0, marker = ".", markersize = 7, capsize=3)
plt.errorbar(3+x, noms(ladungen_4), yerr = stds(ladungen_4), elinewidth = 0.7, linewidth = 0, marker = ".", markersize = 7, capsize=3)
plt.errorbar(4+x, noms(ladungen_5), yerr = stds(ladungen_5), elinewidth = 0.7, linewidth = 0, marker = ".", markersize = 7, capsize=3)

plt.grid()
plt.yticks([0, e, 2*e, 3*e, 4*e, 5*e, 6*e], [0, r"$e$", r"$2e$", r"$3e$", r"$4e$", r"$5e$", r"$6e$"])
plt.xticks([1, 2, 3, 4, 5], [r"$150 \unit{\volt}$", r"$175 \unit{\volt}$", r"$200 \unit{\volt}$", r"$225 \unit{\volt}$", r"$250 \unit{\volt}$"])
plt.ylim(0*e, 6*e)
plt.tight_layout()

plt.savefig("build/plot.pdf")
plt.close()
