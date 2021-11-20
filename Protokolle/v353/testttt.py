import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# # Datenfile erzeugen
# np.random.seed(42)
# N = 50
# data_x = np.linspace(0, 2 * np.pi, N)
# error_y = np.random.lognormal(-1, 0.2, size=N)
# data_y = np.random.normal(np.sin(data_x), error_y)
# np.savetxt(
#     'daten.txt',
#     np.column_stack([data_x, data_y, error_y]),
#     header='x y y_err',
# )

x, y, e_y, d_f = np.genfromtxt('content/data_bc.txt', unpack = True)


def f(x, c):
    return 1/(np.sqrt(1+(x**2 * c**2)))


parameters, pcov = curve_fit(f, x, y)
print(parameters, np.sqrt(np.diag(pcov)), sep='\n')

plt.errorbar(x, y, yerr=e_y, fmt='k.', label='Daten')

t = np.linspace(-0.5, 2 * np.pi + 0.5, 500)
plt.plot(t, f(t, *parameters), label='Fit')
plt.plot(t, np.sin(t), '--', label='Original')

plt.xlim(t[0], t[-1])
plt.xlabel(r'$t$')
plt.ylabel(r'$f(t)$')
plt.legend(loc='best')

plt.tight_layout()
plt.savefig('loesung.pdf')
