#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

head= 'ra,dec,z,c_par,fracDeV_r,velDisp,abs_model_mag_u,abs_model_mag_g,abs_model_mag_r,abs_model_mag_i,abs_model_mag_z,abs_petro_mag_u,abs_petro_mag_g,abs_petro_mag_r,abs_petro_mag_i,abs_petro_mag_z,color_u_r,color_g_r,mu_sup,petroR50_r,petroR90_r'

df = pd.read_csv('/mnt/sda2/extragalactica/2_practico/finals.csv', names=head.split(','), skiprows=1)


print(df.columns)
z= df['z']
M_r= df['abs_model_mag_r']

plt.figure(figsize=(10,6))
plt.scatter( M_r,z, s=1, color='blue')
plt.ylabel('Redshift (z)')
plt.xlabel('Magnitud absoluta (M_r)')
plt.title('Diagrama de Magnitud Absoluta vs Redshift')
plt.plot()

galaxias_brillante = df[df['abs_model_mag_r'] > -15]
print(galaxias_brillante.to_string())

#histogramas de indices de color en un mismo figure
plt.figure(figsize=(10,6))
plt.xlim(0, 3.5)
plt.hist(df['color_u_r'], bins=250, alpha=0.5, label='u-r',edgecolor='blue')
plt.xlabel('Índice de Color (u-r)')
plt.ylabel('Número de Galaxias')
plt.title('Distribución de Índice de Color (u-r)')
plt.legend()
plt.figure(figsize=(10,6))
plt.hist(df['color_g_r'], bins=250, alpha=0.5, label='g-r', color='green')
plt.xlim(-0.6, 1.5)
plt.xlabel('Índice de Color (g-r)')
plt.ylabel('Número de Galaxias')
plt.title('Distribución de Índice de Color (g-r)')
plt.legend()

plt.show()
