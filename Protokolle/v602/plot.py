import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
import scipy.optimize as op


# Konstanten
R = 16.6    #Rydbergenergie in eV
d_lif = 201.4*10**-12 # Gitter in m
alpha = 7.297*10**-3
h = const.h 
e = const.e
c = const.c

# Funktionen
def sigmak(Z, E):                                                   # Funktion zur Berechnung der Abschirmkonstante
    return Z - np.sqrt((E*1000/R) - (((alpha**2)*(Z**4))/4))

def wave_length(theta):
    return 2*d_lif* np.sin(theta*np.pi/180)

def E_K(theta):                                                     # Gibt zu theta zugehörige Energie in eV an
    wave_length = 2*d_lif* np.sin(theta*np.pi/180)
    return h * c / (wave_length*e)


# Ueberprüfen der Bragg-Bedingung
theta_2fach, imps = np.genfromtxt('content/data/Bragg.txt', unpack = True)

theta1 = theta_2fach / 2

plt.plot(theta1, imps, 'x-', color='indianred',label='Messwerte')
plt.plot(27.7/2, 279.0, 'o', color='lightseagreen', label='Maximum', fillstyle='none',ms='6')
plt.axvline(x=27.7/2, color='goldenrod', linestyle = "dotted")
plt.xlabel(r'$\theta \mathbin{/} \unit{\degree}$')
plt.ylabel(r'$\symup{Imp} \mathbin{/} \unit{\second}$')
plt.grid()
plt.legend()

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Bragg.pdf')
plt.close()

print('----------------------')
print(f'Maximum bei Bragg: theta = {27.7/2}° mit 279 imp/s')
print('----------------------')

# Emissionsspektrum
theta_2fach2, imps2 = np.genfromtxt('content/data/Emissionsspektrum.txt', unpack = True)

theta2 = theta_2fach2 / 2

plt.plot(theta2, imps2, '.', color='indianred',label='Messwerte')
plt.plot(theta2[6:80], imps2[6:80], '-', color='dodgerblue',label='Bremsberg')
# plt.plot(10, 800, 'v', color='dodgerblue', label='Bremsberg')
plt.plot(theta2[79:85], imps2[79:85], '-', color='forestgreen', label=r'$K_{\beta}$-Linie')
plt.axvline(x=20.2, color='limegreen', linestyle = "dotted")
# plt.plot(20.2, 1800, 'v', color='forestgreen', label=r'$K_{\beta}$-Linie')
# plt.plot(22.5, 5000, 'v', color='darkorange', label=r'$K_{\alpha}$-Linie')
plt.plot(theta2[89:98], imps2[89:98], '-', color='darkorange', label=r'$K_{\alpha}$-Linie')
plt.axvline(x=22.4, color='burlywood', linestyle = "dotted")
plt.xlabel(r'$\theta \mathbin{/} \unit{\degree}$')
plt.ylabel(r'$\symup{Imp} \mathbin{/} \unit{\second}$')
plt.grid()
plt.legend()

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Emissionsspektrum.pdf')
plt.close()

print('----------------------')
print(f'K_beta Linie:   20.2°')
print(f'K_alpha Linie:  22.4°')
print(f'Beginn Bremsberg: {theta2[6]}°')
print(f'Lambda_min: {wave_length(theta2[6])} m')
print(f'E_min Bremsberg: {E_K(theta2[6])} eV')
print('----------------------')

# Detailspektrum
theta3, imps3 = np.genfromtxt('content/data/Detailspektrum.txt', unpack = True)

plt.plot(theta3, imps3, 'x-', color='indianred',label='Messwerte')
plt.hlines(y=1462/2, xmin=20.05, xmax=20.55, linestyle = ":",color='lightseagreen', linewidth = 1.5,label = r"$\text{Halbwertsbreite des } K_{\beta}\text{-Peaks}$")
plt.hlines(y=4737/2, xmin=22.35, xmax=22.85, linestyle = ":",color='deepskyblue', linewidth = 1.5,label = r"$\text{Halbwertsbreite des } K_{\alpha}\text{-Peaks}$")
plt.xlabel(r'$\theta \mathbin{/} \unit{\degree}$')
plt.ylabel(r'$\symup{Imp} \mathbin{/} \unit{\second}$')
plt.grid()
plt.legend()

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Detailspektrum.pdf')
plt.close()
print('----------------------')
print(f'K_beta Linie Halbwertsbreite:   0.5°')
print(f'K_alpha Linie Halbwertsbreite:  0.5°')
print('----------------------')

# Absorptionsspektren
# Zn
theta_Zn, imps_Zn = np.genfromtxt('content/data/Zn.txt', unpack = True)
plt.plot(theta_Zn, imps_Zn, 'x', color='indianred',label='Messwerte')

# Mittelwerte
unten = np.mean(imps_Zn[0:6])
oben = np.mean(imps_Zn[-5:-1])

plt.hlines(y=unten, xmin= theta_Zn[0], xmax=theta_Zn[6], linestyle=':', linewidth = 1.5, color = 'lightseagreen', label = r'Mittelwert unten')
plt.hlines(y=oben, xmin= theta_Zn[-5], xmax=theta_Zn[-1], linestyle=':', linewidth = 1.5, color = 'deepskyblue', label = r'Mittelwert oben')
plt.axhline(y=(oben+unten)/2, linestyle='--', linewidth = 1.5, color = 'forestgreen', label = r'Mittelwert')
plt.axvline(x=18.6, linestyle='--', linewidth = 1.5, color = 'seagreen', label = r"$\theta_{\symup{K}} = \qty{18,60}{\degree}$")

plt.grid()
plt.legend()
plt.xlabel(r'$\theta \mathbin{/} \unit{\degree}$')
plt.ylabel(r'$\symup{Imp} \mathbin{/} \unit{\second}$')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Zn.pdf')
plt.close()

print('----------------------')
print(f'Theta_K Zn: 18.60°')

# Br
theta_Br, imps_Br = np.genfromtxt('content/data/Br.txt', unpack = True)
plt.plot(theta_Br, imps_Br, 'x', color='indianred',label='Messwerte')

# Mittelwerte
unten = np.mean(imps_Br[0:6])
oben = np.mean(imps_Br[-6:-1])

plt.hlines(y=unten, xmin= theta_Br[0], xmax=theta_Br[6], linestyle=':', linewidth = 1.5, color = 'lightseagreen', label = r'Mittelwert unten')
plt.hlines(y=oben, xmin= theta_Br[-6], xmax=theta_Br[-1], linestyle=':', linewidth = 1.5, color = 'deepskyblue', label = r'Mittelwert oben')
plt.axhline(y=(oben+unten)/2, linestyle='--', linewidth = 1.5, color = 'forestgreen', label = r'Mittelwert')
plt.plot(theta_Br[7:9], imps_Br[7:9], color = "indianred", linewidth = 1, linestyle = "-")
plt.axvline(x=13.18, linestyle='--', linewidth = 1.5, color = 'seagreen', label = r"$\theta_{\symup{K}} = \qty{13,18}{\degree}$")

plt.grid()
plt.legend()
plt.xlabel(r'$\theta \mathbin{/} \unit{\degree}$')
plt.ylabel(r'$\symup{Imp} \mathbin{/} \unit{\second}$')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Br.pdf')
plt.close()

print(f'Theta_K Br: 13.18°')

# Ga
theta_Ga, imps_Ga = np.genfromtxt('content/data/Ga.txt', unpack = True)
plt.plot(theta_Ga, imps_Ga, 'x', color='indianred',label='Messwerte')

# Mittelwerte
unten = np.mean(imps_Ga[0:6])
oben = np.mean(imps_Ga[-6:-1])

plt.hlines(y=unten, xmin= theta_Ga[0], xmax=theta_Ga[6], linestyle=':', linewidth = 1.5, color = 'lightseagreen', label = r'Mittelwert unten')
plt.hlines(y=oben, xmin= theta_Ga[-6], xmax=theta_Ga[-1], linestyle=':', linewidth = 1.5, color = 'deepskyblue', label = r'Mittelwert oben')
plt.axhline(y=(oben+unten)/2, linestyle='--', linewidth = 1.5, color = 'forestgreen', label = r'Mittelwert')
plt.plot(theta_Ga[8:10], imps_Ga[8:10], color = "indianred", linewidth = 1, linestyle = "-")
plt.axvline(x=17.32, linestyle='--', linewidth = 1.5, color = 'seagreen', label = r"$\theta_{\symup{K}} = \qty{17,32}{\degree}$")

plt.grid()
plt.legend()
plt.xlabel(r'$\theta \mathbin{/} \unit{\degree}$')
plt.ylabel(r'$\symup{Imp} \mathbin{/} \unit{\second}$')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Ga.pdf')
plt.close()

print(f'Theta_K Ga: 17.32°')

# Sr
theta_Sr, imps_Sr = np.genfromtxt('content/data/Sr.txt', unpack = True)
plt.plot(theta_Sr, imps_Sr, 'x', color='indianred',label='Messwerte')

# Mittelwerte
unten = np.mean(imps_Sr[0:6])
oben = np.mean(imps_Sr[-7:-1])

plt.hlines(y=unten, xmin= theta_Sr[0], xmax=theta_Sr[6], linestyle=':', linewidth = 1.5, color = 'lightseagreen', label = r'Mittelwert unten')
plt.hlines(y=oben, xmin= theta_Sr[-7], xmax=theta_Sr[-1], linestyle=':', linewidth = 1.5, color = 'deepskyblue', label = r'Mittelwert oben')
plt.axhline(y=(oben+unten)/2, linestyle='--', linewidth = 1.5, color = 'forestgreen', label = r'Mittelwert')
plt.plot(theta_Sr[8:10], imps_Sr[8:10], color = "indianred", linewidth = 1, linestyle = "-")
plt.axvline(x=11.02, linestyle='--', linewidth = 1.5, color = 'seagreen', label = r"$\theta_{\symup{K}} = \qty{11.02}{\degree}$")

plt.grid()
plt.legend()
plt.xlabel(r'$\theta \mathbin{/} \unit{\degree}$')
plt.ylabel(r'$\symup{Imp} \mathbin{/} \unit{\second}$')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Sr.pdf')
plt.close()

print(f'Theta_K Sr: 11.02°')

# Zr
theta_Zr, imps_Zr = np.genfromtxt('content/data/Zr.txt', unpack = True)
plt.plot(theta_Zr, imps_Zr, 'x', color='indianred',label='Messwerte')

# Mittelwerte
unten = np.mean(imps_Zr[0:6])
oben = np.mean(imps_Zr[-6:-1])

plt.hlines(y=unten, xmin= theta_Zr[0], xmax=theta_Zr[6], linestyle=':', linewidth = 1.5, color = 'lightseagreen', label = r'Mittelwert unten')
plt.hlines(y=oben, xmin= theta_Zr[-6], xmax=theta_Zr[-1], linestyle=':', linewidth = 1.5, color = 'deepskyblue', label = r'Mittelwert oben')
plt.axhline(y=(oben+unten)/2, linestyle='--', linewidth = 1.5, color = 'forestgreen', label = r'Mittelwert')
plt.axvline(x=9.90, linestyle='--', linewidth = 1.5, color = 'seagreen', label = r"$\theta_{\symup{K}} = \qty{9,90}{\degree}$")

plt.grid()
plt.legend()
plt.xlabel(r'$\theta \mathbin{/} \unit{\degree}$')
plt.ylabel(r'$\symup{Imp} \mathbin{/} \unit{\second}$')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Zr.pdf')
plt.close()

print(f'Theta_K Zr: 9,90°')
print('----------------------')