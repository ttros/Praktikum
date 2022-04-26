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

plt.plot(theta, imps, '.', color='indianred',label='Messwerte')
plt.plot(theta[0:80], imps[0:80], '-', color='dodgerblue',label='Bremsberg')
# plt.plot(10, 800, 'v', color='dodgerblue', label='Bremsberg')
plt.plot(theta[79:85], imps[79:85], '-', color='forestgreen', label=r'$K_{\beta}$-Linie')
plt.axvline(x=20.2, color='limegreen', linestyle = "dotted")
# plt.plot(20.2, 1800, 'v', color='forestgreen', label=r'$K_{\beta}$-Linie')
# plt.plot(22.5, 5000, 'v', color='darkorange', label=r'$K_{\alpha}$-Linie')
plt.plot(theta[89:98], imps[89:98], '-', color='darkorange', label=r'$K_{\alpha}$-Linie')
plt.axvline(x=22.4, color='burlywood', linestyle = "dotted")
plt.xlabel(r'$\theta \mathbin{/} \unit{\degree}$')
plt.ylabel(r'$\symup{Imp} \mathbin{/} \unit{\second}$')
plt.grid()
plt.legend()

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Emissionsspektrum.pdf')
plt.close()

# Detailspektrum
theta, imps = np.genfromtxt('content/data/Detailspektrum.txt', unpack = True)

plt.plot(theta, imps, 'x-', color='indianred',label='Messwerte')
plt.xlabel(r'$\theta \mathbin{/} \unit{\degree}$')
plt.ylabel(r'$\symup{Imp} \mathbin{/} \unit{\second}$')
plt.grid()
plt.legend()

plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Detailspektrum.pdf')
plt.close()