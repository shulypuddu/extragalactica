#probar repetir esto con el seaborn (https://seaborn.pydata.org/)
# es equivalente a matplotlib pero más sencillo y con mejores gráficos
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv('pr2_shulypuddu.csv')
#print(df.head())

ra = df['ra']
dec = df['dec']
z = df['z']
r_corregido= df['petroMag_r']

plt.figure(figsize=(10, 6))
plt.plot(ra, dec, 'o', markersize=1, color='teal')
plt.xlabel("Ascención recta (ra)")
plt.ylabel("Declinación (dec)")
plt.title("Ubicación de las galaxias en el cielo")
plt.grid(True)
plt.show()


plt.figure(figsize=(10,6))
plt.hist(z, bins=30, color='teal')
plt.xlabel("Redshift (z)")
plt.ylabel("Número de galaxias")
plt.title("Distribución de redshift de las galaxias")
plt.grid(True)
plt.show()


plt.figure(figsize=(10,6))
plt.scatter(z, r_corregido, s=1, color='teal')
plt.ylim(14.5,17.77)
plt.xlabel("Redshift (z)")
plt.ylabel("Magnitud corregida en r (r_corregido)")
plt.title("Magnitud corregida vs Redshift")
plt.grid(True)
plt.show()

