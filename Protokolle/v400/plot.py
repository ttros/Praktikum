import matplotlib.pyplot as plt
import numpy as np

from scipy import stats
from scipy.stats import sem
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
import scipy.constants as const
import scipy.optimize as op

#Aufgabe 2#

alpha_deg = np.array([10,20,30,40,50,60,70])
alpha_rad = ((2*np.pi)/(360))*alpha_deg
beta_deg = np.array([6.5,13,19.5,25.5,31,35.5,39])
beta_rad = ((2*np.pi)/(360))*beta_deg

n_array = np.sin(alpha_rad)/np.sin(beta_rad)
n = ufloat(np.mean(n_array),np.std(n_array))
print(n)
v=2.99793*10**6/n
print(v)