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
#esta galaxia no se xq quedo tan FEA, revisar

#histogramas de indices de color en un mismo figure
fig, axs = plt.subplots(2, 1, figsize=(14, 18)) # 2 filas, 1 columna

# Primer gráfico: scatter de M_r vs z
axs[0].hist(df['color_g_r'], bins=700, alpha=0.7, edgecolor='blue',facecolor='none')
axs[0].set_xlabel('Índice de Color (g-r)')
axs[0].set_xlim(-0.6, 2)
axs[0].set_ylabel('Número de Galaxias')
axs[0].set_title('Histograma color_g_r')

# Segundo gráfico: histograma de color_u_r
axs[1].hist(df['color_u_r'], bins=400, alpha=0.7, edgecolor='green', facecolor='none')
axs[1].set_xlabel('Índice de Color (u-r)')
axs[1].set_xlim(0, 4.5)
axs[1].set_ylabel('Número de Galaxias')
axs[1].set_title('Histograma color_u_r')

plt.tight_layout()
plt.show()

