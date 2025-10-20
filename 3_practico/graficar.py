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
#--------------- Valores preliminares

muestra_10 = galaxias.sample(frac=0.1,random_state=252)
print(f'Tamaño de la muestra: {len(muestra_10)}')
print(galaxias.columns)
fig, axes = plt.subplots(2, 3, figsize=(18, 10))
fig.suptitle('Corrección K a Galaxias ', fontsize=16)

axes[0, 0].axhline(y=0, color='k', linestyle='--', alpha=0.7)
axes[0, 0].scatter(muestra_10['z'], muestra_10['rk_p_u'], marker='.', color='indianred', alpha=0.5, label='Corrección k a $z=0$')
axes[0, 0].scatter(muestra_10['z'], muestra_10['rks_p_u'], marker='.', color='royalblue', alpha=0.5, label='Corrección k a $z=0.1$')
axes[0, 0].axvline(x=0.1, color='k', linestyle='--', alpha=0.7)
axes[0, 0].set_xlabel('$z$')
axes[0, 0].set_ylabel('Corrección K')
axes[0, 0].legend()
axes[0,0].set_title('Filtro $u$')
axes[0, 1].axhline(y=0, color='k', linestyle='--', alpha=0.7)
axes[0, 1].scatter(muestra_10['z'], muestra_10['rk_p_g'], marker='.', color='indianred', alpha=0.5, label='Corrección k a $z=0$')
axes[0, 1].scatter(muestra_10['z'], muestra_10['rks_p_g'], marker='.', color='royalblue', alpha=0.5, label='Corrección k a $z=0.1$')
axes[0, 1].axvline(x=0.1, color='k', linestyle='--', alpha=0.7)
axes[0, 1].set_xlabel('$z$')
axes[0, 1].set_ylabel('Corrección K')
axes[0, 1].legend()
axes[0,1].set_title('Filtro $g$')
axes[0, 2].axhline(y=0, color='k', linestyle='--', alpha=0.7)
axes[0, 2].scatter(muestra_10['z'], muestra_10['rk_p_r'], marker='.', color='indianred', alpha=0.5, label='Corrección k a $z=0$')
axes[0, 2].scatter(muestra_10['z'], muestra_10['rks_p_r'], marker='.', color='royalblue', alpha=0.5, label='Corrección k a $z=0.1$')
axes[0, 2].axvline(x=0.1, color='k', linestyle='--', alpha=0.7)
axes[0, 2].set_xlabel('$z$')
axes[0, 2].set_ylabel('Corrección K')
axes[0, 2].legend()
axes[0, 2].set_title('Filtro $r$')
axes[1, 0].axhline(y=0, color='k', linestyle='--', alpha=0.7)
axes[1, 0].scatter(muestra_10['z'], muestra_10['rk_p_i'], marker='.', color='indianred', alpha=0.5, label='Corrección k a $z=0$')
axes[1, 0].scatter(muestra_10['z'], muestra_10['rks_p_i'], marker='.', color='royalblue', alpha=0.5, label='Corrección k a $z=0.1$')
axes[1, 0].axvline(x=0.1, color='k', linestyle='--', alpha=0.7)
axes[1, 0].set_xlabel('$z$')
axes[1, 0].set_ylabel('Corrección K')
axes[1, 0].legend()
axes[1,0].set_title('Filtro $i$')
axes[1, 1].axhline(y=0, color='k', linestyle='--', alpha=0.7)
axes[1, 1].scatter(muestra_10['z'], muestra_10['rk_p_z'], marker='.', color='indianred', alpha=0.5, label='Corrección k a $z=0$')
axes[1, 1].scatter(muestra_10['z'], muestra_10['rks_p_z'], marker='.', color='royalblue', alpha=0.5, label='Corrección k a $z=0.1$')
axes[1, 1].axvline(x=0.1, color='k', linestyle='--', alpha=0.7)
axes[1, 1].set_xlabel('$z$')
axes[1, 1].set_ylabel('Corrección K')
axes[1, 1].legend()
axes[1,1].set_title('Filtro $z$')
# Remover el subplot vacío
axes[1, 2].remove()
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig('/mnt/sda2/extragalactica/3_practico/informe/imagenes/corrk.pdf', bbox_inches='tight', dpi=300)
#%%
#--------------- Galaxias procesadas

lum_10 = luminosidad.sample(frac=0.1,random_state=252)

print(f"Máximo magnitud absoluta: {luminosidad['petro_abs'].max()}")
print(f"Mínimo magnitud absoluta: {luminosidad['petro_abs'].min()}")
print(f"Volumen máximo: {luminosidad['vmax'].max()}")

plt.axhline(y=0,color='k',linestyle='--', alpha=0.75)
plt.scatter(lum_10['z'],lum_10['vmax'], marker='.',color='royalblue',alpha=0.5)
plt.xlabel('$z$')
plt.ylabel('$V_{max}$ [Mpc$^3$]')
plt.title('Volumen máximo de galaxias en función del corrimiento al rojo')
plt.savefig('/mnt/sda2/extragalactica/3_practico/informe/imagenes/volm.pdf', bbox_inches='tight', dpi=300)
plt.show()
#%%
#--------------- Función de luminosidad 
# Primero normalizo la funcion
M_max=funcion['phi_hist'].max()
phi_hist=funcion['phi_hist']


sns.histplot(luminosidad['petro_abs'],bins='auto',stat='density')
plt.xlabel('Magnitud absoluta')
plt.ylabel(r'$\Phi(M)$')
plt.title('Histograma de magnitudes absolutas')
plt.savefig('/mnt/sda2/extragalactica/3_practico/informe/imagenes/hist_mag_abs.pdf', bbox_inches='tight', dpi=300)
plt.show()
sns.scatterplot(x=funcion['m_medio'],y=phi_hist)
plt.xlabel('Magnitud absoluta')
plt.ylabel(r'$log(\Phi(M))$')
plt.title('Función de luminosidad normalizada')
plt.savefig('/mnt/sda2/extragalactica/3_practico/informe/imagenes/funcion_lum.pdf', bbox_inches='tight', dpi=300)
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

    p0 = [np.log(10.0)*0.4*200, -21.4, -1.2+1]
    popt, cov = opt.curve_fit(schechter, funcion['m_medio'], funcion['phi_hist']*np.diff(funcion['m_medio'])[0], p0=p0)
    print(f'{popt=}')

    return popt, cov


popt, cov = fit_lumfunc()
# Graficar los datos y la función ajustada
M_fit = np.linspace(funcion['m_medio'].min(), funcion['m_medio'].max())
phi_fit = schechter(M_fit, *popt)+12
plt.plot(M_fit, phi_fit, label='Ajuste Schechter', color='red')
sns.scatterplot(x=funcion['m_medio'],y=phi_hist)
plt.xlabel('Magnitud absoluta')
plt.ylabel(r'$\Phi(M)$')
plt.title('Función de luminosidad normalizada')
plt.savefig('/mnt/sda2/extragalactica/3_practico/informe/imagenes/funcion_lum.pdf', bbox_inches='tight', dpi=300)
plt.show()








