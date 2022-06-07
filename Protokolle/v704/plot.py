import matplotlib.pyplot as plt
import numpy as np
import uncertainties as unc
from uncertainties import unumpy as unp 
from uncertainties import ufloat
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as devs
import scipy.optimize as op


def linfit(x,m,b):
    return m*x+b

def regression(x,y,x1,x2,Element,Farbe,Name):
    params, pcov = op.curve_fit(linfit, x, y)
    std = np.sqrt(np.diag(pcov))
    m = params[0]
    b = params[1]
    M=unp.uarray(m,std[0])
    B=unp.uarray(b,std[1])

    print(f'Lambda {Element}: {-M}')
    print(f'N_0 {Element}: {unp.exp(B)}')
    print()

    xx = np.linspace(x1,x2,1000)
    plt.plot(xx,m*xx+b,color = Farbe, label=Name)
    return M,B

# background_per_sec = 184/600

# ####Vanadium####

# counts_vn_raw = np.genfromtxt('content/data/data_Vn.txt', unpack = True)
# counts_vn=counts_vn_raw-(background_per_sec*30)
# delta_vn = np.around(np.sqrt(counts_vn), 0)

# log_counts_vn = np.log(counts_vn)
# t_vn=np.linspace(30,900,30)

# plt.errorbar(t_vn, log_counts_vn, xerr=None, yerr=delta_vn/counts_vn, fmt='.', color='tomato',label='Messwerte')
# m_vn,b_vn = regression(t_vn,log_counts_vn,0,900,'Vanadium','dodgerblue','Regression')
# plt.xlabel(r'$t \mathbin{/} \unit{\second}$')
# plt.ylabel(r'$\symup{log}\,N$')
# plt.grid()
# plt.legend()
# plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
# plt.savefig('build/Vanadium.pdf')
# plt.close()
# print(f'Halbwertszeit Vanadium: {np.log(2)/-m_vn}')


# #####Rhodium######

# print()
# print()

# counts_rh_raw = np.genfromtxt('content/data/data_Rh.txt', unpack = True)
# counts_rh=counts_rh_raw-(background_per_sec*15)
# delta_rh = np.around(np.sqrt(counts_rh), 0)

# log_counts_rh = np.log(counts_rh)
# t_rh=np.linspace(15,720,48)

# plt.errorbar(t_rh, log_counts_rh, xerr=None, yerr=delta_rh/counts_rh, fmt='.', color='tomato',label='Messwerte')
# t_long=t_rh[20:]
# plt.vlines(t_long[0], 2, 7, 'darkslategray', 'dashed',label='$t^*$')
# array_long = log_counts_rh[20:]
# m_rh_long,b_rh_long = regression(t_long, array_long,t_long[0],720,'Rhodium 104','dodgerblue','Regression für $\ce{^104 Rh}$')

# t_short = t_rh[:15]
# counts_rh_short=counts_rh[:15]
# for i in range(15):
#     counts_rh_short[i]=counts_rh_short[i]-np.exp(noms(m_rh_long)*t_rh[i]+noms(b_rh_long))

# log_counts_rh_short = np.log(counts_rh_short)
# m_rh_short,b_rh_short = regression(t_short, log_counts_rh_short,0,t_short[14],'Rhodium 104i','gold','Regression für $\ce{^{104\symup{i}} Rh}$')

# Summenkurve = np.exp(noms(m_rh_short)*t_rh+noms(b_rh_short)) + np.exp(noms(m_rh_long)*t_rh+noms(b_rh_long))
# plt.plot(t_rh, np.log(Summenkurve), '-', color='indigo', label='berechnete Summenkurve')

# plt.xlabel(r'$t \mathbin{/} \unit{\second}$')
# plt.ylabel(r'$\symup{log}\,N$')
# plt.grid()
# plt.legend()
# plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
# plt.savefig('build/Rhodium.pdf')
# plt.close()
# print(f'Halbwertszeit Rhodium 104: {np.log(2)/-m_rh_long}')
# print(f'Halbwertszeit Rhodium 104i: {np.log(2)/-m_rh_short}')