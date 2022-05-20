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
gamma=60
n_kron=1.510
alpha_1=np.array([30.0,35.0,40.0,45.0,50.0])
alpha_2_rot=np.array([79.0,69.0,60.0,53.0,47.0])
alpha_2_grün=np.array([81.0,70.0,61.0,54.0,48.0])

beta_1_rad=np.arcsin((np.sin((2*np.pi)/(360)*alpha_1))/n_kron)
beta_1=beta_1_rad*((360)/(2*np.pi))
beta_2=gamma-beta_1

delta_rot = (alpha_1+alpha_2_rot)-(beta_1+beta_2)
delta_grün = (alpha_1+alpha_2_grün)-(beta_1+beta_2)

print(f'Ablenkung rot: {delta_rot}')
print(f'Ablenkung grün: {delta_grün}')

#Aufgabe 5#
def lambda_(d,k,phi):
    lambda_k = d*np.sin(((2*np.pi)/(360))*phi)/k
    return lambda_k

print(f'600mm rot: {lambda_((1/600)*10**6,1,24)}')
print(f'600mm grün: {lambda_((1/600)*10**6,1,19.5)}')

print()

print(f'300mm rot k=1: {lambda_((1/300)*10**6,1,11)}')
print(f'300mm rot k=2: {lambda_((1/300)*10**6,2,23)}')
print(f'300mm grün k=1: {lambda_((1/300)*10**6,1,9.5)}')
print(f'300mm grün k=2: {lambda_((1/300)*10**6,2,19)}')
print(f'300mm grün k=3: {lambda_((1/300)*10**6,3,29)}')

print()

print(f'100mm rot k=1: {lambda_((1/100)*10**6,1,4)}')
print(f'100mm rot k=2: {lambda_((1/100)*10**6,2,6.5)}')
print(f'100mm rot k=3: {lambda_((1/100)*10**6,3,11)}')
print(f'100mm rot k=4: {lambda_((1/100)*10**6,4,15)}')
print(f'100mm rot k=5: {lambda_((1/100)*10**6,5,19.5)}')
print(f'100mm grün k=1: {lambda_((1/100)*10**6,1,4)}')
print(f'100mm grün k=2: {lambda_((1/100)*10**6,2,6.5)}')
print(f'100mm grün k=3: {lambda_((1/100)*10**6,3,9.5)}')
print(f'100mm grün k=4: {lambda_((1/100)*10**6,4,12.5)}')
print(f'100mm grün k=5: {lambda_((1/100)*10**6,5,16)}')

print()

lambda_r = np.array([lambda_((1/600)*10**6,1,24),lambda_((1/300)*10**6,1,11),lambda_((1/300)*10**6,2,23),lambda_((1/100)*10**6,1,4),lambda_((1/100)*10**6,2,6.5),lambda_((1/100)*10**6,3,11),lambda_((1/100)*10**6,4,15),lambda_((1/100)*10**6,5,19.5)])
lambda_rot  = ufloat(np.mean(lambda_r),np.std(lambda_r))
print(f'rot mittel: {lambda_rot}')

lambda_g = np.array([lambda_((1/600)*10**6,1,19.5),lambda_((1/300)*10**6,1,9.5),lambda_((1/300)*10**6,2,19),lambda_((1/300)*10**6,3,29),lambda_((1/100)*10**6,1,4),lambda_((1/100)*10**6,2,6.5),lambda_((1/100)*10**6,3,9.5),lambda_((1/100)*10**6,4,12.5),lambda_((1/100)*10**6,5,16)])
lambda_gruen  = ufloat(np.mean(lambda_g),np.std(lambda_g))
print(f'grün mittel: {lambda_gruen}')