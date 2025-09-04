
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy 
import astropy.units as u
from astropy.cosmology import Planck18 as cosmo


#cosmo.H0 = 70 * cosmo.H0.unit
#cosmo.Om0 = 0.3 # Matter density parameter
#cosmo.Ode0 = 0.7 # Dark energy density parameter

comov_dist=[]
dist_lum =[]
dist_ang =[]
z = np.arange(0,5,0.001)

for i in range(len(z)):
    comov_dist = cosmo.comoving_distance(i*0.001)
    dist_lum = cosmo.luminosity_distance(i*0.001)
    dist_ang = cosmo.angular_diameter_distance(i*0.001)

plt.figure(figsize=(10, 6))
plt.plot(z, comov_coord)
plt.xlabel("Redshift")
plt.ylabel("Comoving Distance (Mpc)")
plt.title("Comoving Distance vs. Redshift")
plt.grid(True)


plt.figure(figsize=(10, 6))
plt.plot(z, dist_lum)
plt.xlabel("Redshift")
plt.ylabel("Distancia luminosa (Mpc)")
plt.title("Distancia luminosa vs. Redshift")
plt.grid(True)


plt.figure(figsize=(10, 6))
plt.plot(z, dist_ang)
plt.xlabel("Redshift")
plt.ylabel("Distancia angular (Mpc)")
plt.title("Distancia angular vs. Redshift")
plt.grid(True)
plt.show()

df = pd.read_csv('cosmos.csv')
print(df.head())


