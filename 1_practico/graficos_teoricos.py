#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import astropy.units as u
from astropy.cosmology import WMAP5 as cosmo

t= np.arange(0,1,0.0001)
H0=cosmo.H0.value
def universo_polvoso(t):
    a = (t**(2/3))*2/(3*H0)
    return a

def universo_radiativo(t):
    a = (t**(1/2))*2*H0
    return a
def universo_lamb(t):
    a= np.exp(H0*t)   
    return a


a_mat = universo_polvoso(t)
a_rad = universo_radiativo(t)
a_lamb = universo_lamb(t)

# Grafica de factor de escala de materia
plt.figure(figsize=(10, 6))
plt.plot(t,a_mat, color='orange', label='Materia')
plt.plot(t, a_rad, color='teal', label='Radiación')
plt.plot(t, a_lamb, color='green', label='Materia Oscura')
plt.xlabel("Tiempo t")
plt.ylabel("Factor de escala a(t)")
plt.title("Evolución del factor de escala en el universo")
plt.ylim(0,10)


plt.show()
