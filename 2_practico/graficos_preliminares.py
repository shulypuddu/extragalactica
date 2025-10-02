#probar repetir esto con el seaborn (https://seaborn.pydata.org/)
# es equivalente a matplotlib pero más sencillo y con mejores gráficos
#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv('/mnt/sda2/extragalactica/2_practico/pr2_shulypuddu.csv')

#print(df.head())
bar_color="#7166D1"
edge_color= "#695ED1"
linea_color='tomato'

ra = df['ra']
dec = df['dec']
z = df['z']
r_corregido= df['petroMag_r']-df['extinction_r']

plt.figure(figsize=(10, 6))
plt.plot(ra, dec, 'o', markersize=1, color=bar_color)
plt.xlabel("Ascención recta (ra)")
plt.ylabel("Declinación (dec)")
plt.title("Ubicación de las galaxias en el cielo")
plt.grid(True)
plt.savefig('/mnt/sda2/extragalactica/2_practico/imagenes/gal_in_sky.png', dpi=300,bbox_inches='tight')

plt.figure(figsize=(10,6))
plt.hist(z, bins=30, edgecolor=edge_color,facecolor=bar_color)
plt.xlabel("Redshift (z)")
plt.ylabel("Número de galaxias")
plt.title("Distribución de redshift de las galaxias")
plt.savefig('/mnt/sda2/extragalactica/2_practico/imagenes/gal_z.png', dpi=300,bbox_inches='tight')


plt.figure(figsize=(10,6))
plt.scatter(z, r_corregido, s=1, color=bar_color)
plt.ylim(14.5,17.77)
plt.xlabel("Redshift (z)")
plt.ylabel("Magnitud corregida en r (r_corregido)")
plt.title("Magnitud corregida vs Redshift")
plt.savefig('/mnt/sda2/extragalactica/2_practico/imagenes/mag_vs_z.png', dpi=300,bbox_inches='tight')

# %%
