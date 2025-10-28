#%%
#--------------- Librerias del programa y cargar los datos.
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy.stats as st
import scipy.optimize as opt


galaxias = pd.read_fwf('/mnt/sda2/extragalactica/3_practico/datos.dat')
head_gal = 'z,petro_r,red_r,r50,rk_p_u,rk_p_g,rk_p_r,rk_p_i,rk_p_z,rks_p_u,rks_p_g,rks_p_r,rks_p_i,rks_p_z'
galaxias.columns= head_gal.split(',')

luminosidad = pd.read_fwf('/mnt/sda2/extragalactica/3_practico/valores.dat',skiprows=1,names=[ 'z', 'dL', 'petro_abs', 'rks_p', 'dLmax', 'zmax', 'vmax' ])
funcion = pd.read_fwf('/mnt/sda2/extragalactica/3_practico/fun_lum.dat',skiprows=1,names=[ 'm_medio', 'phi', 'phi_hist' ])
print(luminosidad.columns)
print(funcion.columns)

#%% 
#--------------- Corrección K a Galaxias

muestra_10 = galaxias.sample(frac=0.1,random_state=252)
print(f'Tamaño de la muestra: {len(muestra_10)}')
print(galaxias.columns)


# Crear una sola figura con 2 filas y 5 columnas
fig, axes = plt.subplots(2, 5, figsize=(22, 8))
fig.suptitle('Corrección K a Galaxias', fontsize=16)

# Primera fila - Corrección K a z=0
for i in range(5):
    filtros = ['u', 'g', 'r', 'i', 'z']
    columnas = ['rk_p_u', 'rk_p_g', 'rk_p_r', 'rk_p_i', 'rk_p_z']
    
    axes[0, i].axhline(y=0, color='k', linestyle='--', alpha=0.7)
    axes[0, i].scatter(muestra_10['z'], muestra_10[columnas[i]], marker='.', color='indianred', alpha=0.5)
    axes[0, i].axvline(x=0.0, color='k', linestyle='--', alpha=0.7)
    axes[0, i].set_xlabel('$z$')
    axes[0, i].set_title(f'Filtro ${filtros[i]}$')
    if i == 0:
        axes[0, i].set_ylabel('Corrección K (z=0)')

# Segunda fila - Corrección K a z=0.1
for i in range(5):
    filtros = ['u', 'g', 'r', 'i', 'z']
    columnas = ['rks_p_u', 'rks_p_g', 'rks_p_r', 'rks_p_i', 'rks_p_z']
    
    axes[1, i].axhline(y=0, color='k', linestyle='--', alpha=0.7)
    axes[1, i].scatter(muestra_10['z'], muestra_10[columnas[i]], marker='.', color='royalblue', alpha=0.5)
    axes[1, i].axvline(x=0.1, color='k', linestyle='--', alpha=0.7)
    axes[1, i].set_xlabel('$z$')
    axes[1, i].set_title(f'Filtro ${filtros[i]}$')
    if i == 0:
        axes[1, i].set_ylabel('Corrección K (z=0.1)')

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig('/mnt/sda2/extragalactica/3_practico/informe/imagenes/corrk_combined.pdf', bbox_inches='tight', dpi=300)
plt.show()

#%%
# Crear una sola figura con 5 filas y 1 columna
fig, axes = plt.subplots(5, 1, figsize=(22, 40))
fig.suptitle('Corrección K a Galaxias a z=0', fontsize=16)

# Primera fila - Corrección K a z=0
for i in range(5):
    filtros = ['u', 'g', 'r', 'i', 'z']
    columnas = ['rk_p_u', 'rk_p_g', 'rk_p_r', 'rk_p_i', 'rk_p_z']
    
    axes[i].axhline(y=0, color='k', linestyle='--', alpha=0.7)
    axes[i].scatter(muestra_10['z'], muestra_10[columnas[i]], marker='.', color='indianred', alpha=0.5)
    axes[i].axvline(x=0.0, color='k', linestyle='--', alpha=0.7)
    axes[i].set_xlabel('$z$')
    axes[i].set_title(f'Filtro ${filtros[i]}$')

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig('/mnt/sda2/extragalactica/3_practico/informe/imagenes/corrk_0.png', bbox_inches='tight', dpi=300)
plt.show()


fig, axes = plt.subplots(5, 1, figsize=(22, 40))
fig.suptitle('Corrección K a Galaxias a z=0.1', fontsize=16)

for i in range(5):
    filtros = ['u', 'g', 'r', 'i', 'z']
    columnas = ['rks_p_u', 'rks_p_g', 'rks_p_r', 'rks_p_i', 'rks_p_z']
    
    axes[i].axhline(y=0, color='k', linestyle='--', alpha=0.7)
    axes[i].scatter(muestra_10['z'], muestra_10[columnas[i]], marker='.', color='royalblue', alpha=0.5)
    axes[i].axvline(x=0.1, color='k', linestyle='--', alpha=0.7)
    axes[i].set_xlabel('$z$')
    axes[i].set_title(f'Filtro ${filtros[i]}$')

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig('/mnt/sda2/extragalactica/3_practico/informe/imagenes/corrk_1.png', bbox_inches='tight', dpi=300)
plt.show()








#%%
#--------------- Galaxias procesadas

lum_10 = luminosidad.sample(frac=0.1,random_state=252)

print(f"Máximo magnitud absoluta: {luminosidad['petro_abs'].max()}")
print(f"Mínimo magnitud absoluta: {luminosidad['petro_abs'].min()}")
print(f"Volumen máximo: {luminosidad['vmax'].max()}")


plt.scatter(lum_10['z'],lum_10['vmax'], marker='.',color='royalblue',alpha=0.5)
plt.xlabel('$z$')
plt.ylabel('$V_{max}$ [Mpc$^3$]')
plt.title('Volumen máximo de galaxias en función del corrimiento al rojo')
plt.savefig('/mnt/sda2/extragalactica/3_practico/informe/imagenes/volm.pdf', bbox_inches='tight', dpi=300)
plt.show()

sns.histplot(x=galaxias['z'],bins=300,color='indianred',alpha=0.7,stat='density')
plt.xlim(-0.005,0.37)
plt.xlabel('$z$')
plt.ylabel('Número de galaxias')
plt.title('Histograma de corrimiento al rojo')
plt.axvline(x=0.1, color='k', linestyle='--', alpha=0.7)
plt.savefig('/mnt/sda2/extragalactica/3_practico/informe/imagenes/hist_z.pdf', bbox_inches='tight', dpi=300)
plt.show()

#%%
#--------------- Función de luminosidad 
# Primero normalizo la funcion
M_max=funcion['phi_hist'].max()
phi_hist=funcion['phi_hist']


sns.histplot(luminosidad['petro_abs'],bins=20,stat='density',color='royalblue',alpha=0.7)
plt.xlabel('Magnitud absoluta')
plt.ylabel(r'$\Phi(M)$')
plt.title('Histograma de magnitudes absolutas')
plt.savefig('/mnt/sda2/extragalactica/3_practico/informe/imagenes/hist_mag_abs.pdf', bbox_inches='tight', dpi=300)
plt.show()
sns.scatterplot(x=funcion['m_medio'],y=phi_hist,color='royalblue')
plt.xlabel('Magnitud absoluta')
plt.ylabel(r'$log(\Phi(M))$')
plt.title('Función de luminosidad normalizada')
plt.savefig('/mnt/sda2/extragalactica/3_practico/informe/imagenes/histo_lum.pdf', bbox_inches='tight', dpi=300)
plt.show()

#%%
def schechter(x, logA, B, C):
    """
    fit func para log(phi)
    x = M: mag absoulta
    logA = phi^{*} ln(10) : ...
    B = M^{*} : magnitud de escala
    C = 1+alpha : pendiente
    returns log10(phi) 
    """
    return logA - 0.4*(x-B)*C - 10.0**(-0.4*(x-B))/np.log(10.0)

def fit_lumfunc():

    p0 = [np.log(10.0)*0.4*M_max, -21.4, -1.2+1]
    popt, cov = opt.curve_fit(schechter, funcion['m_medio'], funcion['phi_hist'], p0=p0)
    print(f'{popt=}')

    return popt, cov


popt, cov = fit_lumfunc()
# Graficar los datos y la función ajustada
M_fit = np.linspace(funcion['m_medio'].min(), funcion['m_medio'].max())
phi_fit = schechter(M_fit, *popt)
plt.plot(M_fit, phi_fit, label='Ajuste Schechter', color='indianred')
sns.scatterplot(x=funcion['m_medio'],y= funcion['phi_hist'],color='royalblue', label='Función observada')
plt.xlabel('Magnitud absoluta')
plt.ylabel(r'$\Phi(M)$')
plt.title('Función de luminosidad normalizada')
plt.savefig('/mnt/sda2/extragalactica/3_practico/informe/imagenes/funcion_lum.pdf', bbox_inches='tight', dpi=300)
plt.show()

#%%





