import matplotlib.pyplot as plt
import numpy as np

import uncertainties as unp
from scipy import stats
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
import scipy.constants as const

c = 2.99792458*10**8

def regression(U,I,x1,x2,name):
    m , b , r ,p ,std =stats.linregress(U,I)
    M=unp.uarray(m,std)
    B=unp.uarray(b,std)

    print(f'm {name}: {M}')
    print(f'b {name}: {B}')
    print()

    xx = np.linspace(x1,x2,1000)
    plt.plot(xx,m*xx+b,color = "dodgerblue", label='Fit')
    return m,b

def Ugrenz(m,b,farbe):
    U=-b/m
    print(f'Grenzspannung {farbe}: {U}')
    print()
    return U

U_r, I_r = np.genfromtxt('content/data/rot.txt', unpack = True)
plt.plot(U_r, np.sqrt(I_r),'x', color='indianred',label='Messwerte')
plt.xlabel(r'$U_{\symup{B}}\,/\,\unit{\volt}$')
plt.ylabel(r'$\sqrt{I}\,/\,\sqrt{\unit{\ampere}}$')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
m_r, b_r = regression(U_r,np.sqrt(I_r),-0.4,0.7,'rot')
plt.grid()
plt.legend(loc='best')
plt.savefig('build/plot_rot.pdf')
plt.close()

U_g, I_g = np.genfromtxt('content/data/gelb.txt', unpack = True)
plt.plot(U_g, np.sqrt(I_g),'x', color='indianred',label='Messwerte')
plt.xlabel(r'$U_{\symup{B}}\,/\,\unit{\volt}$')
plt.ylabel(r'$\sqrt{I}\,/\,\sqrt{\unit{\ampere}}$')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
m_g, b_g = regression(U_g,np.sqrt(I_g),0.2,0.5,'gelb')
plt.grid()
plt.legend(loc='best')
plt.savefig('build/plot_gelb.pdf')
plt.close()

U_gr, I_gr = np.genfromtxt('content/data/gruen.txt', unpack = True)
plt.plot(U_gr, np.sqrt(I_gr),'x', color='indianred',label='Messwerte')
plt.xlabel(r'$U_{\symup{B}}\,/\,\unit{\volt}$')
plt.ylabel(r'$\sqrt{I}\,/\,\sqrt{\unit{\ampere}}$')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
m_gr, b_gr = regression(U_gr,np.sqrt(I_gr),0.1,0.6,'gruen')
plt.grid()
plt.legend(loc='best')
plt.savefig('build/plot_gruen.pdf')
plt.close()

U_v1, I_v1 = np.genfromtxt('content/data/violett_1.txt', unpack = True)
plt.plot(U_v1, np.sqrt(I_v1),'x', color='indianred',label='Messwerte')
plt.xlabel(r'$U_{\symup{B}}\,/\,\unit{\volt}$')
plt.ylabel(r'$\sqrt{I}\,/\,\sqrt{\unit{\ampere}}$')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
m_v1, b_v1 = regression(U_v1,np.sqrt(I_v1),0.0,1.2,'violett 1')
plt.grid()
plt.legend(loc='best')
plt.savefig('build/plot_violett_1.pdf')
plt.close()

U_v2, I_v2 = np.genfromtxt('content/data/violett_2.txt', unpack = True)
plt.plot(U_v2, np.sqrt(I_v2),'x', color='indianred',label='Messwerte')
plt.xlabel(r'$U_{\symup{B}}\,/\,\unit{\volt}$')
plt.ylabel(r'$\sqrt{I}\,/\,\sqrt{\unit{\ampere}}$')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
m_v2, b_v2 = regression(U_v2,np.sqrt(I_v2),0.0,1.3,'violett 2')
plt.grid()
plt.legend(loc='best')
plt.savefig('build/plot_violett_2.pdf')
plt.close()

U_g_rot = Ugrenz(m_r,b_r, 'rot')
U_g_gelb = Ugrenz(m_g,b_g,'gelb')
U_g_gruen = Ugrenz(m_gr,b_gr,'gruen')
U_g_violett1 = Ugrenz(m_v1,b_v2,'violett 1')
U_g_violett2 = Ugrenz(m_v2,b_v2,'violett 2')

#U_g = np.array([U_g_rot,U_g_gelb,U_g_gruen,U_g_violett1,U_g_violett2])
#lambda_=np.array([646*10**(-9),578*10**(-9),546*10**(-9),435*10**(-9),405*10**(-9)])
#f=c/lambda_

##wegen kagg daten##############################################################################
U_g_fix = np.array([U_g_gelb,U_g_gruen,U_g_violett1,U_g_violett2])
lambda_fix=np.array([578*10**(-9),546*10**(-9),435*10**(-9),405*10**(-9)])
################################################################################################
f_fix=c/lambda_fix

plt.plot(f_fix, U_g_fix, 'x', color='indianred',label='Grenzspannungen')
plt.plot(c/(646*10**(-9)),U_g_rot, color='dimgray', marker='x')
plt.xlabel(r'$f\,/\,\unit{\hertz}$')
plt.ylabel(r'$U_G\,/\,\unit{\volt}$')
mm, bb = regression(f_fix,U_g_fix,0,8*10**14,'Fit')
plt.grid()
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_planck.pdf')
plt.close()

U_gg, I_gg = np.genfromtxt('content/data/gelb_lang.txt', unpack = True)
plt.plot(U_gg,I_gg,'x-', color='indianred',label='Messwerte')
plt.xlabel(r'$U\,/\,\unit{\volt}$')
plt.ylabel(r'$I\,/\,\unit{\ampere}$')
plt.grid()
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_gelb_lang.pdf')
plt.close()