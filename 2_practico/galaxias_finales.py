#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as st
import scipy.optimize as opt

#nombre de las columnas
head= 'ra,dec,z,c_par,fracDeV_r,velDisp,abs_model_mag_u,abs_model_mag_g,abs_model_mag_r,abs_model_mag_i,abs_model_mag_z,abs_petro_mag_u,abs_petro_mag_g,abs_petro_mag_r,abs_petro_mag_i,abs_petro_mag_z,color_u_r,color_g_r,mu_sup,petroR50_r,petroR90_r'

#leo el archivo. Paso el head como los nombres aclarando que el nombre de cada columna esta separado por una coma.
df0 = pd.read_csv('/mnt/sda2/extragalactica/2_practico/finals.csv', names=head.split(','), skiprows=1)

galaxias_brillante = df0[df0['abs_model_mag_r'] > -15]
#esta galaxia no se xq quedo tan FEA, asique la filtro
df = df0[df0['abs_model_mag_r'] < -15]


print(df.columns)
z= df['z']
M_r= df['abs_model_mag_r']
#%%
plt.figure(figsize=(10,6))
plt.scatter( M_r,z, s=1, color='teal', marker='.')
plt.ylabel('Redshift ($z$)')
plt.xlabel('Magnitud absoluta ($M_r$)')
plt.title('Diagrama de Magnitud Absoluta vs Redshift')
plt.plot()


#%% ----------histogramas de indices de color en un mismo figure
fig, axs = plt.subplots(2, 1, figsize=(14, 18)) # 2 filas, 1 columna

# Primer gráfico: scatter de M_r vs z
axs[0].hist(df['color_g_r'], bins='auto', alpha=0.55,edgecolor='teal',facecolor='teal')
axs[0].set_xlabel('Índice de Color (g-r)')
axs[0].set_xlim(-0.6, 2)
axs[0].set_ylabel('Número de Galaxias')
axs[0].set_title('Histograma color_g_r')

# Segundo gráfico: histograma de color_u_r
axs[1].hist(df['color_u_r'], bins='auto', alpha=0.55,edgecolor='hotpink', facecolor='hotpink')
axs[1].set_xlabel('Índice de Color (u-r)')
axs[1].set_xlim(0, 4.5)
axs[1].set_ylabel('Número de Galaxias')
axs[1].set_title('Histograma color_u_r')

plt.tight_layout()
plt.show()

#Ahora voy a ajustar una doble gaussiana a cada figura
u_r=st.norm.fit(df['color_u_r'])
g_r=st.norm.fit(df['color_g_r'])
fig, axs = plt.subplots(2, 1, figsize=(14, 18)) # 2 filas, 1 columna

# Primer gráfico: scatter de M_r vs z
axs[0].hist(df['color_g_r'], bins='auto', alpha=0.55,edgecolor='teal',facecolor='teal')
axs[0].set_xlabel('Índice de Color (g-r)')
axs[0].set_xlim(-0.6, 2)
axs[0].set_ylabel('Número de Galaxias')
axs[0].set_title('Histograma color_g_r')
axs[0].plot(np.linspace(-0.6,2,100), st.norm.pdf(np.linspace(-0.6,2,100), g_r[0], g_r[1])*len(df['color_g_r'])*0.5, color='darkblue', label='Ajuste Gaussiano')
axs[0].legend()
# Segundo gráfico: histograma de color_u_r
axs[1].hist(df['color_u_r'], bins='auto', alpha=0.55,edgecolor='hotpink', facecolor='hotpink')
axs[1].set_xlabel('Índice de Color (u-r)')
axs[1].set_xlim(0, 4.5)
axs[1].set_ylabel('Número de Galaxias')
axs[1].set_title('Histograma color_u_r')

plt.tight_layout()
plt.show()



#%%
#-----------------------------------------------------------------
# Histograma de las distribuciones del parámetro de concentración y de fracDeV
plt.figure(figsize=(10,6))
plt.hist(df['c_par'], bins='auto',edgecolor='hotpink',facecolor='none', density=True)
plt.xlim(0,6)
plt.title('Distribución del parámetro de concentración (C)')

plt.figure(figsize=(10,6))
plt.hist(df['fracDeV_r'], bins='auto',edgecolor='teal',facecolor='none', density=True)
plt.xlim(-0.0001,1.0001)
plt.title('Distribución del parámetro fracDeV en el filtro r')

#%% Creo que este es el que pide en el inciso donde habla de correlación
#Todo esto hay que hacerlo ya separando segun si son de la nube azul o de la secuencia roja 

plt.figure(figsize=(10,6))
plt.scatter(df['c_par'],df['color_u_r'], marker='.')

# DIAGRAMA COLOR MAGNITUD
plt.figure(figsize=(10,6))
plt.scatter(df['color_u_r'], df['abs_petro_mag_r'],marker='.',color='hotpink')
plt.xlim(-1,13)

#%%

plt.figure(figsize=(10,6))
plt.scatter(df['abs_model_mag_r'], df['petroR50_r'], s=1, color='teal', marker='.')
plt.ylabel('Radio que contiene el 50% de la luz en r (petroR50_r)')
plt.xlabel('Magnitud absoluta en r ($M_r$)')
plt.title('Diagrama de Magnitud Absoluta vs Radio (50%)')
plt.plot()

