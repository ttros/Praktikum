import numpy as np
import uncertainties as unp
from scipy import stats
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)

def Mittelwert(array,n):
    Mittel=array.mean()
    Standartabweichung=array.std(ddof=1)
    Mittelwert=ufloat(Mittel,Standartabweichung/np.sqrt(n))
    return Mittelwert

n=10
g=9.81



#### Kurze Länge ####

l=0.5

T_pendel1_5x, T_pendel2_5x, T_gleichphasig_5x, T_gegenphasig_5x, T_gekoppelt_5x, T_schwebung= np.genfromtxt('content/data/data_kurz.txt', unpack = True)

T_gekoppelt_5x=T_gekoppelt_5x[:-1]
T_schwebung=T_schwebung[:-1]

T_pendel1 = T_pendel1_5x/5
T_pendel2 = T_pendel2_5x/5
T_gleichphasig = T_gleichphasig_5x/5
T_gegenphasig = T_gegenphasig_5x/5
T_gekoppelt = T_gekoppelt_5x/5

print(f'----------Kurzes Pendel----------')
print('\n')
print(f'Mittelwert 1. Pendel: {Mittelwert(T_pendel1,n)} s')
print(f'Mittelwert 2. Pendel: {Mittelwert(T_pendel2,n)} s')
print('\n')

print(f'Mittelwert gleichphasig: {Mittelwert(T_gleichphasig,n)} s')
T_gleichphasig_theorie = 2*np.pi*np.sqrt(l/g)
print(f'gleichphasig theoretisch: {T_gleichphasig_theorie} s')
print('\n')

print(f'Mittelwert gegenphasig: {Mittelwert(T_gegenphasig,n)} s')
kappa=( Mittelwert(T_gleichphasig,n)**2 - Mittelwert(T_gegenphasig,n)**2 )/( Mittelwert(T_gleichphasig,n)**2 + Mittelwert(T_gegenphasig,n)**2 )
print(f'Kappa: {kappa}')
T_gegenphasig_theorie = 2*np.pi*unp.sqrt(l/(g+2*kappa))
print(f'gegenphasig theoretisch: {T_gegenphasig_theorie} s')
print('\n')

print(f'Mittelwert gekoppelt: {Mittelwert(T_gekoppelt,n-1)} s')
print(f'Mittelwert Schwebung: {Mittelwert(T_schwebung,n-1)} s')
T_schwebung_theorie = ( Mittelwert(T_gleichphasig,n) * Mittelwert(T_gegenphasig,n) )/( Mittelwert(T_gleichphasig,n) - Mittelwert(T_gegenphasig,n) )
print(f'Schwebung theoretisch: {T_schwebung_theorie} s')
print('\n')



#### Lange Länge ####

l=1

T_pendel1_5x, T_pendel2_5x, T_gleichphasig_5x, T_gegenphasig_5x, T_gekoppelt_5x, T_schwebung= np.genfromtxt('content/data/data_lang.txt', unpack = True)

T_pendel1 = T_pendel1_5x/5
T_pendel2 = T_pendel2_5x/5
T_gleichphasig = T_gleichphasig_5x/5
T_gegenphasig = T_gegenphasig_5x/5
T_gekoppelt = T_gekoppelt_5x/5

print(f'----------Langes Pendel----------')
print('\n')
print(f'Mittelwert 1. Pendel: {Mittelwert(T_pendel1,n)} s')
print(f'Mittelwert 2. Pendel: {Mittelwert(T_pendel2,n)} s')
print('\n')

print(f'Mittelwert gleichphasig: {Mittelwert(T_gleichphasig,n)} s')
T_gleichphasig_theorie = 2*np.pi*np.sqrt(l/g)
print(f'gleichphasig theoretisch: {T_gleichphasig_theorie} s')
print('\n')

print(f'Mittelwert gegenphasig: {Mittelwert(T_gegenphasig,n)} s')
kappa=( Mittelwert(T_gleichphasig,n)**2 - Mittelwert(T_gegenphasig,n)**2 )/( Mittelwert(T_gleichphasig,n)**2 + Mittelwert(T_gegenphasig,n)**2 )
print(f'Kappa: {kappa}')
T_gegenphasig_theorie = 2*np.pi*unp.sqrt(l/(g+2*kappa))
print(f'gegenphasig theoretisch: {T_gegenphasig_theorie} s')
print('\n')

print(f'Mittelwert gekoppelt: {Mittelwert(T_gekoppelt,n)} s')
print(f'Mittelwert Schwebung: {Mittelwert(T_schwebung,n)} s')
T_schwebung_theorie = ( Mittelwert(T_gleichphasig,n) * Mittelwert(T_gegenphasig,n) )/( Mittelwert(T_gleichphasig,n) - Mittelwert(T_gegenphasig,n) )
print(f'Schwebung theoretisch: {T_schwebung_theorie} s')
print('\n')