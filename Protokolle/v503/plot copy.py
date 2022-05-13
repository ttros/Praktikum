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
d = ufloat(7.5250*10**-3,0.0051*10**-3)     # in m
dichte_oel = 886                            # in kg/m^3
dichte_L = 1.204                            # in kg/m^3 bei 20°C

# %%% Definition Formeln %%%

def ladung(v_auf, v_ab, U, eta_L):
    return 3*np.pi*eta_L*np.sqrt((9*eta_L*(v_ab-v_auf))/4*const.g*(dichte_oel-dichte_L))*(v_ab+v_auf)/(U/d)

def radius(v_auf, v_ab, eta_L):
    return np.sqrt((9*eta_L*(v_ab-v_auf))/2*const.g*(dichte_oel-dichte_L))


# times = ['t_auf_1', 't_ab_1', 't_auf_2', 't_ab_2', 't_auf_3', 't_ab_3', 't_auf_4', 't_ab_4', 't_auf_5', 't_ab_5']
j=1
t_auf_1, t_ab_1, t_auf_2, t_ab_2, t_auf_3, t_ab_3, t_auf_4, t_ab_4, t_auf_5, t_ab_5 = np.genfromtxt(f'content/data/Messung_{j}.txt', unpack = True)

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
print(times('t_ab_1'))
# velocities_mean = {
#     'v_m_auf_1': ufloat(np.mean(velocities('v_auf_1')),sem(velocities('v_auf_1'))),
#     'v_m_ab_1' : ufloat(np.mean(velocities('v_ab_1' )),sem(velocities('v_ab_1' ))),
#     'v_m_auf_2': ufloat(np.mean(velocities('v_auf_2')),sem(velocities('v_auf_2'))),
#     'v_m_ab_2' : ufloat(np.mean(velocities('v_ab_2' )),sem(velocities('v_ab_2' ))),
#     'v_m_auf_3': ufloat(np.mean(velocities('v_auf_3')),sem(velocities('v_auf_3'))),
#     'v_m_ab_3' : ufloat(np.mean(velocities('v_ab_3' )),sem(velocities('v_ab_3' ))),
#     'v_m_auf_4': ufloat(np.mean(velocities('v_auf_4')),sem(velocities('v_auf_4'))),
#     'v_m_ab_4' : ufloat(np.mean(velocities('v_ab_4' )),sem(velocities('v_ab_4' ))),
#     'v_m_auf_5': ufloat(np.mean(velocities('v_auf_5')),sem(velocities('v_auf_5'))),
#     'v_m_ab_5' : ufloat(np.mean(velocities('v_ab_5' )),sem(velocities('v_ab_5' )))
# }

for x in times:
    print(times[f'{x}'])

for x in velocities:
    print(velocities[f'{x}'])