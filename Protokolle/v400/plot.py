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

alpha_deg = np.array([10.0,20.0,30.0,40.0,50.0,60.0,70.0])
alpha_rad = ((2*np.pi)/(360))*alpha_deg
beta_deg = np.array([6.5,13.0,19.5,25.5,31.0,35.5,39.0])
beta_rad = ((2*np.pi)/(360))*beta_deg

n_array = np.sin(alpha_rad)/np.sin(beta_rad)
n = ufloat(np.mean(n_array),np.std(n_array))
print(f'Brechungsindex: {n}')
v=2.99793*10**6/n
print(f'c in Plexiglas: {v}')

#Aufgabe 3#

d=5.85
s_methode1=d*(np.sin(alpha_rad-beta_rad))/(np.cos(beta_rad))
print(f'Strahlenversatz 1.Methode: {s_methode1}')
beta_calc=np.arcsin(np.sin(alpha_rad)/np.mean(n_array))
s_methode2=d*(np.sin(alpha_rad-beta_calc))/(np.cos(beta_calc))
print(f'Strahlenversatz 2.Methode: {s_methode2}')

#Aufgabe 4#
gamma=((2*np.pi)/(360))*60
n_kron=1.510
alpha_1=np.array([30.0,35.0,40.0,45.0,50.0])
alpha_2_rot=np.array([79.0,69.0,60.0,53.0,47.0])
alpha_2_grün=np.array([81.0,70.0,61.0,54.0,48.0])

beta_1=np.arcsin(np.sin(alpha_1)/n_kron)
beta_2=gamma-beta_1

delta_rot = (alpha_1+alpha_2_rot)-(beta_1+beta_2)
delta_grün = (alpha_1+alpha_2_grün)-(beta_1+beta_2)

print(f'Ablenkung rot: {delta_rot}')
print(f'Ablenkung grün: {delta_grün}')