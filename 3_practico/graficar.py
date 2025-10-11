#%%
#--------------- Librerias del programa y cargar los datos.
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

galaxias = pd.read_table('/mnt/sda2/extragalactica/3_practico/datos.dat', sep=r'\s+')
head_gal = 'z,petro_r,red_r,r50,rk_p_u,rk_p_g,rk_p_r,rk_p_i,rk_p_z,rks_p_u,rks_p_g,rks_p_r,rks_p_i,rks_p_z'
galaxias.columns= head_gal.split(',')

luminosidad = pd.read_csv('/mnt/sda2/extragalactica/3_practico/fort.csv')
funcion = pd.read_csv('/mnt/sda2/extragalactica/3_practico/fun_lum.csv')
funcion.columns

#%% 
#--------------- Valores preliminares

muestra_10 = galaxias.sample(frac=0.1,random_state=252)
print(f'Tamaño de la muestra: {len(muestra_10)}')
print(galaxias.columns)
plt.axhline(y=0,color='k')
plt.scatter(muestra_10['z'],muestra_10['rk_p_r'], marker='.',color='indianred',alpha=0.5)
plt.scatter(muestra_10['z'],muestra_10['rks_p_r'], marker='.',color='royalblue',alpha=0.5)

#%%
#--------------- Galaxias procesadas

lum_10 = luminosidad.sample(frac=0.1,random_state=252)

print(f"Máximo magnitud absoluta: {luminosidad['petro_abs'].max()}")
print(f"Mínimo magnitud absoluta: {luminosidad['petro_abs'].min()}")
print(f"Volumen máximo: {luminosidad['vmax'].max()}")

plt.axhline(y=0,color='k')
plt.scatter(lum_10['z'],lum_10['vmax'], marker='.',color='b',alpha=0.5)

#%%
#--------------- Función de luminosidad 
sns.histplot(luminosidad['petro_abs'],bins='auto')
plt.show()
sns.scatterplot(x=funcion['m_medio'],y=funcion['phi_hist'])


