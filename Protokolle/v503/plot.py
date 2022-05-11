import matplotlib.pyplot as plt
import numpy as np

import uncertainties as unp
from scipy import stats
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
