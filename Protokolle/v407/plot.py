import matplotlib.pyplot as plt
import numpy as np
import uncertainties as unc
from uncertainties import unumpy as unp 
from uncertainties import ufloat
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as stds
import scipy.optimize as op
import scipy.constants as const

def n_s_pol(alpha, E):
    return np.sqrt((2*E*np.cos(2*alpha)+ E**2 + 1)/(1-2*E + E**2))

def n_p_pol(alpha, E):
    b = ((E+1)/(E-1))**2
    return np.sqrt(b/(2*np.cos(alpha)**2) + np.sqrt(b**2/(4*np.cos(alpha)**4) - b*np.tan(alpha)**2))

def fresnel_s(alpha, n):
    return -(np.cos(alpha) - np.sqrt(n**2-np.sin(alpha)**2)) / (np.cos(alpha) + np.sqrt(n**2-np.sin(alpha)**2))

def fresnel_p(alpha, n):
    return (n**2*np.cos(alpha) - np.sqrt(n**2-np.sin(alpha)**2)) / (n**2*np.cos(alpha) + np.sqrt(n**2-np.sin(alpha)**2)) 

alpha, I_s_raw, I_p_raw = np.genfromtxt('content/data/data.txt', unpack = True)

I_0_s=(540/3)*10**-6
I_0_p=(520/3)*10**-6
I_dunkel=7*10**-9

I_s_raw*=10**-6
I_p_raw*=10**-6

I_s = I_s_raw - I_dunkel
I_p = I_p_raw - I_dunkel

n_s = n_s_pol((alpha*2*np.pi)/360, np.sqrt((I_s)/I_0_s))
n_p = n_p_pol((alpha*2*np.pi)/360, np.sqrt((I_p)/I_0_p))

n_s_mean = np.mean(n_s[n_s > 3])
n_s_std = np.std(n_s[n_s > 3])
n_s_mittelwert = ufloat(n_s_mean, n_s_std)
print(f'Mittelwert für n s-polarisation: {n_s_mittelwert}')

n_p_mean = np.mean(n_p[n_p < 4.4])
n_p_std = np.std(n_p[n_p < 4.4])
n_p_mittelwert = ufloat(n_p_mean, n_p_std)
print(f'Mittelwert für n p-polarisation: {n_p_mittelwert}')

plt.plot(alpha, np.sqrt((I_s)/I_0_s),'x', color='tomato', label='Messwerte s-polarisiertes Licht')
plt.plot(alpha, np.sqrt((I_p)/I_0_p),'x', color='steelblue', label='Messwerte p-polarisiertes Licht')
plt.vlines(75, ymin = 0, ymax = 1, color='darkolivegreen', linestyle = ':', label='Brewsterwinkel')

alpha_lin = np.linspace(5,87,1000)
alpha_lin_rad = ((alpha_lin*2*np.pi)/360)
plt.plot(alpha_lin,fresnel_s(alpha_lin_rad, n_s_mean), color='green', label='Theoriekurve s-polarisiertes Licht')
brewster=73.5
alpha_lin_1=np.linspace(5,brewster,1000)
alpha_lin_1_rad = ((alpha_lin_1*2*np.pi)/360)
alpha_lin_2=np.linspace(brewster,87,500)
alpha_lin_2_rad = ((alpha_lin_2*2*np.pi)/360)
plt.plot(alpha_lin_1,fresnel_p(alpha_lin_1_rad, n_p_mean), color='gold', label='Theoriekurve p-polarisiertes Licht')
plt.plot(alpha_lin_2,-fresnel_p(alpha_lin_2_rad, n_p_mean), color='gold')

plt.xlabel(r'$\alpha \mathbin{/} \unit{\degree}$')
plt.ylabel(r'$\sqrt{\frac{I(\alpha)}{I_0}}$')
plt.grid()
plt.legend()
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot.pdf')
plt.close()

brewster=75
n = np.tan(brewster*2*np.pi/360)
print(f'Brechungsindex: {n}')