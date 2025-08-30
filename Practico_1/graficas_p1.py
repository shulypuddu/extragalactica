#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy 
import astropy.units as u
from astropy.cosmology import FlatLambaCDM as cosmo


cosmo.H0 = 70 * cosmo.H0.unit



z = np.arange(0,5,0.001)
comov_dist = cosmo.comoving_distance(z)
dist_lum = cosmo.luminosity_distance(z)
dist_ang = cosmo.angular_diameter_distance(z)

plt.figure(figsize=(10, 6))
plt.plot(z, comoving_coord)
plt.xlabel("Redshift")
plt.ylabel("Comoving Distance (Mpc)")
plt.title("Comoving Distance vs. Redshift")
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(z, dist_lum)
plt.xlabel("Redshift")
plt.ylabel("Distancia luminosa (Mpc)")
plt.title("Distancia luminosa vs. Redshift")
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(z, dist_ang)
plt.xlabel("Redshift")
plt.ylabel("Distancia angular (Mpc)")
plt.title("Distancia angular vs. Redshift")
plt.grid(True)
plt.show()

df = pd.read_csv('cosmos.csv')
print(df.head())
