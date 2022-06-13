import matplotlib.pyplot as plt
import numpy as np
import uncertainties as unc
from uncertainties import unumpy as unp 
from uncertainties import ufloat
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as stds
import scipy.optimize as op
import scipy.constants as const


def linfit(x,m,b):
    return m*x+b

def regression(x,y,x1,x2,Farbe,Name):
    params, pcov = op.curve_fit(linfit, x, y)
    std = np.sqrt(np.diag(pcov))
    m = params[0]
    b = params[1]
    M=unp.uarray(m,std[0])
    B=unp.uarray(b,std[1])

    xx = np.linspace(x1,x2,1000)
    plt.plot(xx,m*xx+b,color = Farbe, label=Name)
    return M,B

