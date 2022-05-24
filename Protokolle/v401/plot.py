import matplotlib.pyplot as plt
import numpy as np
import uncertainties as unc
from uncertainties import unumpy as unp 
from uncertainties import ufloat
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as devs


n1 = np.genfromtxt("content/data/data_1.txt", unpack = True)
n2 = np.genfromtxt("content/data/data_2.txt", unpack = True)

u = 1/5.017         # Untersetzungsverhältnis
d = 5e-3            # m Verschiebung des Spiegels 
laser = 635e-9      # m Wellenlänge des Lasers
b = 50e-3           # m Schichtdicke der Messzelle

T_0 = 273.15        # K
T   = 296.14        # K, Raumtemperatur (23°C)
p_0 = 101325        # Pascal
d_p =  20000        # Pascal Druch nach Vakkumpumpen
b = 0.05            # Dicke der Druckkammer in m


z1 = ufloat(np.mean(n1), np.std(n1))
print("Mittelwert der Zählraten: ", z1)
print("Armlängenänderung", '{0:.5e}'.format(d*u))
laser_exp = 2*d*u/z1
d_laser = np.abs(noms(laser_exp)- laser)/laser

print("Experimentelle Wellenlänge: ", '{0:.5e}'.format(laser_exp))
print("Abweichung:                 ", d_laser)

z2 = ufloat(np.mean(n2), np.std(n2))
print("Mittelwert der Zählraten(Druck): ", z2)
print("Delta(n):                        ", z2*laser/(2*b))
n_exp = 1 + z2*laser/(2*b)*(T/T_0)*(p_0/(p_0 - d_p))
print("Experimenteller Brechungsindex: ", '{0:.6f}'.format(n_exp))

print(n1)