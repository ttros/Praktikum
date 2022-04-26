import matplotlib.pyplot as plt
import numpy as np

# Ueberprüfen der Bragg-Bedingung
theta_2fach, imps = np.genfromtxt('content/data/Bragg.txt', unpack = True)

theta = theta_2fach / 2

plt.plot(theta, imps, 'x', color='indianred',label='Messwerte')
plt.plot(27.7/2, 279.0, 'o', color='lightseagreen', label='Maximum', fillstyle='none',ms='6')
plt.axvline(x=27.7/2, color='goldenrod', linestyle = "dotted")
plt.xlabel(r'$\theta \mathbin{/} \unit{\degree}$')
plt.ylabel(r'$\symup{Imp} \mathbin{/} \unit{\second}$')
plt.grid()
plt.legend()

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Bragg.pdf')
plt.close()

# Emissionsspektrum
theta_2fach, imps = np.genfromtxt('content/data/Emissionsspektrum.txt', unpack = True)

theta = theta_2fach / 2

plt.plot(theta, imps, 'x-', color='indianred',label='Messwerte')
plt.plot(10, 8000, 'v', color='dodgerblue', label='Bremsberg')
plt.xlabel(r'$\theta \mathbin{/} \unit{\degree}$')
plt.ylabel(r'$\symup{Imp} \mathbin{/} \unit{\second}$')
plt.grid()
plt.legend()

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Emissionsspektrum.pdf')
