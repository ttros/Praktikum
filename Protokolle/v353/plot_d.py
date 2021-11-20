import matplotlib.pyplot as plt
import numpy as np
import uncertainties as unp
from scipy.optimize import curve_fit

f, A, a, b = np.genfromtxt('content/data_bc.txt', unpack = True)
phi = a/b * 2 * np.pi

def f1(f,c):
    return np.arctan(-f*c)

parameters, pcov = curve_fit(f1, f, phi, sigma=None)
print(parameters, np.sqrt(np.diag(pcov)), sep='\n')

phi2 = np.linspace(0, (np.pi)/2, 50)

plt.polar(phi, A, 'rx')
plt.polar(phi2, np.sin(phi2)/(np.tan(phi2))*2.8, 'b')

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_d.pdf')
