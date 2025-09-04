#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import astropy.units as u
from astropy.cosmology import Planck18 as cosmo

t= np.arange(0,1e4,1)
H0= 70e-6 #(70 km/s/Mpc) en 1/Myr
def universo_polvoso(t):
    a = (t*2/(3)*H0)**(2/3)
    return a

def universo_radiativo(t):
    a = (t*2*H0)**(1/2)
    return a
def universo_lamb(t):
    a= np.exp(H0*t)   
    return a

a_mat = universo_polvoso(t)
a_rad = universo_radiativo(t)
a_lamb = universo_lamb(t)

# Grafica de factor de escala de materia
plt.figure(figsize=(10, 6))
plt.plot(t,a_mat, color='purple')
plt.xlabel("Tiempo [yr]")
plt.ylabel("Factor de escala de materia")
plt.title("materia vs. Tiempo")
plt.grid(True)

# Grafica de factor de escala de radiacion
plt.figure(figsize=(10, 6))
plt.plot(t, a_rad, color='purple')
plt.xlabel("Tiempo [yr]")
plt.ylabel("Factor de escala de radiacion")
plt.title("radiacion vs. tiempo")
plt.grid(True)

# Grafica de factor de escala de energia oscura
plt.figure(figsize=(10, 6))
plt.plot(t, a_lamb, color='purple')
plt.xlabel("Tiempo [yr]")
plt.ylabel("Factor de escala de energia oscura")
plt.title("energia oscura vs. tiempo")
plt.grid(True)


plt.show()

# %%
