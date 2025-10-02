#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as st
import scipy.optimize as opt

#nombre de las columnas
head= 'ra,dec,z,c_par,fracDeV_r,velDisp,abs_model_mag_u,abs_model_mag_g,abs_model_mag_r,abs_model_mag_i,abs_model_mag_z,abs_petro_mag_u,abs_petro_mag_g,abs_petro_mag_r,abs_petro_mag_i,abs_petro_mag_z,color_u_r,color_g_r,mu_sup,petroR50_r,petroR90_r'

#leo el archivo. Paso el head como los nombres aclarando que el nombre de cada columna esta separado por una coma.
df0 = pd.read_csv('/mnt/sda2/extragalactica/2_practico/galaxias_finales.csv', names=head.split(','), skiprows=1)
bar_color="#7166D1"
edge_color= "#695ED1"
linea_color='tomato'

galaxias_brillante = df0[df0['abs_model_mag_r'] > -15]
#esta galaxia no se xq quedo tan FEA, asique la filtro
df = df0[df0['abs_model_mag_r'] < -15]

color_u_r = df['color_u_r'].dropna()

# Definir la función suma de dos gaussianas
def doble_gaussiana(x, amp1, mu1, sigma1, amp2, mu2, sigma2):
    return (amp1 * np.exp(-(x - mu1)**2 / (2 * sigma1**2)) +
            amp2 * np.exp(-(x - mu2)**2 / (2 * sigma2**2)))

def gauss(x,amp,mu,sigma):
    return amp * np.exp(-(x - mu)**2 / (2 * sigma**2))

# Obtener histograma (cuentas y bordes)

cuentas, bordes = np.histogram(color_u_r, bins=1000,density=True)
centros = (bordes[:-1] + bordes[1:]) / 2

# Estimaciones iniciales para los parámetros
p0 = [1, color_u_r.mean()-0.5, color_u_r.std()/2, 1, color_u_r.mean()+0.5, color_u_r.std()/2]

# Ajuste
params, cov = opt.curve_fit(doble_gaussiana, centros, cuentas, p0=p0)

x_fit = np.linspace(color_u_r.min(), color_u_r.max(), 200)
x1_fit = np.linspace(color_u_r.min(), 3, 200)
x2_fit = np.linspace( 1,color_u_r.max(), 200)
g_azul= gauss(x1_fit,params[0],params[1],params[2])
g_rojo= gauss(x2_fit,params[3],params[4],params[5])

def diferencia(x):
    return gauss(x,params[0],params[1],params[2]) - gauss(x,params[3],params[4],params[5])
x_inicial= (params[4]+params[1])/2
x_cruce = opt.fsolve(diferencia, x_inicial)[0]

rojo = df[df['color_u_r'] > x_cruce]
azul = df[df['color_u_r'] < x_cruce]

tempranas = df[df['c_par'] >= 2.5]
tardias = df[df['c_par'] < 2.5]

bulge = df[df['fracDeV_r'] >= 0.8]
disco = df[df['fracDeV_r']<= 0.2]

def separar_cuartiles_numpy(datos):
    """
    Separa datos en cuartiles usando NumPy
    """
    datos = np.array(datos)
    
    # Calcular los valores que dividen en cuartiles
    q1 = np.percentile(datos, 25)
    q2 = np.percentile(datos, 50)  # Mediana
    q3 = np.percentile(datos, 75)
    
    # Separar los datos en los 4 cuartiles
    cuartil_1 = datos[datos <= q1]
    cuartil_2 = datos[(datos > q1) & (datos <= q2)]
    cuartil_3 = datos[(datos > q2) & (datos <= q3)]
    cuartil_4 = datos[datos > q3]
    
    return {
        'valores_division': [q1, q2, q3],
        'cuartiles': {
            'Q1': cuartil_1,
            'Q2': cuartil_2, 
            'Q3': cuartil_3,
            'Q4': cuartil_4
        },
        'tamaños': [len(cuartil_1), len(cuartil_2), len(cuartil_3), len(cuartil_4)]
    }
mag_r_quartiles=separar_cuartiles_numpy(df['abs_petro_mag_r'])

def ajuste_lineal(x,y):
    """
    --Parámetros--
    x, y : tipo array - lista de datos para ajustar
    
    
    --Retorna-- 
    Lista de parametros y errores

    """
    _y=np.log(y)
    param,pcov = np.polyfit(x,_y,1,cov=True)

    return param[0]*x+param[1]


lineal_rojo= ajuste_lineal(rojo['abs_model_mag_r'],rojo['petroR50_r'])
lineal_azul= ajuste_lineal(azul['abs_model_mag_r'],azul['petroR50_r'])
lineal_tempranas= ajuste_lineal(tempranas['abs_model_mag_r'],tempranas['petroR50_r'])
lineal_tardias= ajuste_lineal(tardias['abs_model_mag_r'],tardias['petroR50_r'])
lineal_bulge= ajuste_lineal(bulge['abs_model_mag_r'],bulge['petroR50_r'])
lineal_disco= ajuste_lineal(disco['abs_model_mag_r'],disco['petroR50_r'])



#%% ---------Diagrama de redshift vs magnitud absoluta
plt.figure(figsize=(10,6))
plt.scatter(df['z'],df['abs_model_mag_r'], s=1, color=bar_color, marker='.',alpha=0.75)
plt.xlabel('Redshift ($z$)', fontsize=14)
plt.ylabel('Magnitud absoluta ($M_r$)', fontsize=14)
plt.title('Diagrama de Redshift vs Magnitud Absoluta',fontsize=17)
plt.tight_layout()
plt.gca().invert_yaxis()  
plt.savefig('/mnt/sda2/extragalactica/2_practico/imagenes/mabs_vs_z.png', dpi=300,bbox_inches='tight')

plt.show()  

#%% ----------histogramas de indices de color en un mismo figure
fig, axs = plt.subplots(2, 1, figsize=(14, 18)) # 2 filas, 1 columna

# Primer gráfico: scatter de M_r vs z
axs[0].hist(df['color_g_r'], bins='auto', edgecolor=edge_color,facecolor=bar_color,density=True)
axs[0].set_xlabel('Índice de Color (g-r)',fontsize='large')
axs[0].set_xlim(-0.25, 1.25)
axs[0].set_ylabel('Número de Galaxias',fontsize='large')
axs[0].set_title('Histograma color_g_r',fontsize=16)

# Segundo gráfico: histograma de color_u_r
axs[1].hist(df['color_u_r'], bins='auto',edgecolor=edge_color, facecolor=bar_color,density=True)
axs[1].set_xlabel('Índice de Color (u-r)',fontsize='large')
axs[1].set_xlim(0, 3.65)
axs[1].set_ylabel('Número de Galaxias',fontsize='large')
axs[1].set_title('Histograma color_u_r',fontsize=16)

plt.tight_layout()

plt.savefig('/mnt/sda2/extragalactica/2_practico/imagenes/colores.png', dpi=300,bbox_inches='tight')

plt.show()

#Ahora voy a ajustar una doble gaussiana a cada figura
#%%
# Supón que tus datos están en la variable color_u_r

# Graficar
plt.hist(color_u_r, bins='auto', density=True, edgecolor=edge_color, facecolor=bar_color, label='Datos')
plt.plot(x_fit, doble_gaussiana(x_fit, *params), label='Ajuste Doble Gaussiana', color='darkslategray')
plt.plot(x1_fit, g_azul,  label='Doble Gaussiana',color='dodgerblue') 
plt.plot(x2_fit, g_rojo, color='indianred', label='Doble Gaussiana') 
plt.xlim(0, 4)
plt.xlabel('color_u_r', fontsize=14) 
plt.ylabel('Densidad', fontsize=14)
plt.title('Ajuste gaussiano a cada pico de color',fontsize=17)
plt.savefig('/mnt/sda2/extragalactica/2_practico/imagenes/colores_ajustados.png', dpi=300,bbox_inches='tight')

plt.show()

#%% -------------Ahora separo por galaxias rojas y azules

fig, axs = plt.subplots(2, 1, figsize=(14, 18)) # 2 filas, 1 columna
# Primer gráfico: scatter de M_r vs z
axs[0].hist(rojo['color_g_r'], bins='auto',density=True, alpha=0.55,facecolor='indianred',edgecolor='indianred')
axs[0].hist(azul['color_g_r'], bins='auto',density=True, alpha=0.55,facecolor='dodgerblue',edgecolor='dodgerblue')
axs[0].set_xlabel('Índice de Color (g-r)', fontsize=14)
axs[0].set_xlim(0, 1)
axs[0].set_ylabel('Número de Galaxias', fontsize=14)
axs[0].set_title('Histograma color_g_r', fontsize=17)
# Segundo gráfico: histograma de color_u_r
axs[1].hist(rojo['color_u_r'], bins='auto',density=True, alpha=0.55,facecolor='indianred',edgecolor='indianred')
axs[1].hist(azul['color_u_r'], bins='auto',density=True, alpha=0.55,facecolor='dodgerblue',edgecolor='dodgerblue')
axs[1].set_xlabel('Índice de Color (u-r)', fontsize=14)
axs[1].set_xlim(0.5,3.25)
axs[1].set_ylabel('Número de Galaxias', fontsize=14)
axs[1].set_title('Histograma color_u_r', fontsize=17)

plt.savefig('/mnt/sda2/extragalactica/2_practico/imagenes/rojas_y_azules.png', dpi=300,bbox_inches='tight')
plt.tight_layout()
plt.show()

#%% ----Para definir los cuartiles (a simple vista) hago un histograma con 4 bines para la abs. mag.
plt.hist(df['abs_model_mag_r'],bins=4,density=True)


#%%
#-----------------------------------------------------------------
# Histograma de las distribuciones del parámetro de concentración y de fracDeV
plt.figure(figsize=(10,6))
plt.hist(df['c_par'], bins='auto',edgecolor=edge_color,facecolor=bar_color, density=True)
plt.xlim(1.15,4.0)
plt.title('Distribución del parámetro de concentración (C)')
plt.savefig('/mnt/sda2/extragalactica/2_practico/imagenes/hist_c.png', dpi=300,bbox_inches='tight')

#fracDeV habla de cuanto bulge hay en la galaxia.

plt.figure(figsize=(10,6))
plt.hist(df['fracDeV_r'], bins='auto',edgecolor=edge_color,facecolor=bar_color, density=True)
plt.xlim(-0.0001,1.0001)
plt.title('Distribución del parámetro fracDeV en el filtro r')
plt.savefig('/mnt/sda2/extragalactica/2_practico/imagenes/hist_fracdev.png', dpi=300,bbox_inches='tight')


plt.figure(figsize=(10,6))
plt.scatter(df['c_par'],df['fracDeV_r'], marker='.',color=bar_color,alpha=0.5)
plt.xlabel('Parámetro de concentración (C)')
plt.ylabel('Frac de V')
plt.xlim(0.75,4)
plt.savefig('/mnt/sda2/extragalactica/2_practico/imagenes/fracdev_vs_c.png', dpi=300,bbox_inches='tight')

#%% Creo que este es el que pide en el inciso donde habla de correlación
#Todo esto hay que hacerlo ya separando segun si son de la nube azul o de la secuencia roja 

#plt.title()

plt.figure(figsize=(10,6))
plt.scatter(df['color_u_r'],df['c_par'], marker='.',color=bar_color)
#plt.scatter(rojo['color_u_r'],rojo['c_par'],marker='.',color='indianred',alpha=0.5)
#plt.scatter(azul['color_u_r'],azul['c_par'],marker='.',color='dodgerblue',alpha=0.5)
plt.ylabel('(C)')
plt.xlabel('Indice de color u-r')
plt.xlim(0.25,3.5)
plt.ylim(1.25,4.15)
plt.savefig('/mnt/sda2/extragalactica/2_practico/imagenes/ur_vs_c.png', dpi=300,bbox_inches='tight')


plt.figure(figsize=(10,6))
plt.scatter(df['color_g_r'],df['c_par'], marker='.',color=bar_color)
plt.scatter(rojo['color_g_r'],rojo['c_par'],marker='.',color='indianred',alpha=0.5)
plt.scatter(azul['color_g_r'],azul['c_par'],marker='.',color='dodgerblue',alpha=0.5)
plt.ylabel('(C)')
plt.xlabel('Indice de color g-r')
plt.xlim(-0.25,1.5)
plt.ylim(1.25,4.15)
plt.savefig('/mnt/sda2/extragalactica/2_practico/imagenes/gr_vs_c.png', dpi=300,bbox_inches='tight')


#%% DIAGRAMA COLOR MAGNITUD
plt.figure(figsize=(10,6))
plt.scatter(df['color_u_r'],df['abs_petro_mag_r'],  marker='.',color=bar_color,alpha=0.5)
plt.scatter(rojo['color_u_r'],rojo['abs_petro_mag_r'],  marker='.',color='indianred',alpha=0.5)
plt.scatter(azul['color_u_r'],azul['abs_petro_mag_r'], marker='.',color='dodgerblue',alpha=0.5)
plt.xlim(0,4)
plt.xlabel('Índice de color (u-r)')
plt.ylabel('$M_r$')
plt.title('Diagrama color - magnitud')
plt.gca().invert_yaxis()  # Invertir eje Y (magnitud)
plt.plot()
plt.savefig('/mnt/sda2/extragalactica/2_practico/imagenes/dcm_ur.png', dpi=300,bbox_inches='tight')
plt.show()

plt.figure(figsize=(10,6))
plt.scatter(df['color_g_r'],df['abs_petro_mag_r'],  marker='.',color=bar_color,alpha=0.5)
plt.scatter(rojo['color_g_r'],rojo['abs_petro_mag_r'],  marker='.',color='indianred',alpha=0.5)
plt.scatter(azul['color_g_r'],azul['abs_petro_mag_r'], marker='.',color='dodgerblue',alpha=0.5)
plt.xlim(-0.25,1.25)
plt.xlabel('Índice de color (g-r)')
plt.ylabel('$M_r$')
plt.title('Diagrama color - magnitud')
plt.gca().invert_yaxis()  # Invertir eje Y (magnitud)
plt.plot()
plt.savefig('/mnt/sda2/extragalactica/2_practico/imagenes/dcm_gr.png', dpi=300,bbox_inches='tight')
plt.show()



#%%

plt.figure(figsize=(10,6))
plt.scatter(df['abs_model_mag_r'], df['petroR50_r'], s=1, color=bar_color, marker='.')
plt.yscale('log')
plt.ylim(1,10)
plt.ylabel('Radio petrosiano 50$\%$')
plt.xlabel('Magnitud absoluta en r ($M_r$)')
plt.title('Diagrama de Magnitud Absoluta vs Radio ($50\%$)')
plt.savefig('/mnt/sda2/extragalactica/2_practico/imagenes/mabs_vs_r50.png', dpi=300,bbox_inches='tight')

plt.figure(figsize=(10,6))
plt.scatter(rojo['abs_model_mag_r'], rojo['petroR50_r'], s=1, color='indianred', marker='.',alpha=0.5)
plt.scatter(azul['abs_model_mag_r'], azul['petroR50_r'], s=1, color='dodgerblue', marker='.',alpha=0.5)
plt.yscale('log')
plt.plot(azul['abs_model_mag_r'],lineal_azul,color='dodgerblue')
plt.ylim(1,10)
plt.ylabel('Radio petrosiano 50$\%$')
plt.xlabel('Magnitud absoluta en r ($M_r$)')
plt.title('Diagrama de Magnitud Absoluta vs Radio ($50\%$)')
plt.savefig('/mnt/sda2/extragalactica/2_practico/imagenes/mabs_vs_r50_azulyrojo.png', dpi=300,bbox_inches='tight')

plt.figure(figsize=(10,6))
plt.scatter(tempranas['abs_model_mag_r'],tempranas['petroR50_r'], s=1, color='teal', marker='.',alpha=0.5)
plt.scatter(tardias['abs_model_mag_r'], tardias['petroR50_r'], s=1, color='hotpink', marker='.',alpha=0.5)
plt.yscale('log')
plt.ylim(1,10)
plt.ylabel('Radio petrosiano 50$\%$')
plt.xlabel('Magnitud absoluta en r ($M_r$)')
plt.title('Diagrama de Magnitud Absoluta vs Radio ($50\%$)')
plt.savefig('/mnt/sda2/extragalactica/2_practico/imagenes/mabs_vs_r50_earlyylate.png', dpi=300,bbox_inches='tight')


plt.figure(figsize=(10,6))
plt.scatter(bulge['abs_model_mag_r'],bulge['petroR50_r'], s=1, color='darkorange', marker='.',alpha=0.5)
plt.scatter(disco['abs_model_mag_r'], disco['petroR50_r'], s=1, color='royalblue', marker='.',alpha=0.5)
plt.yscale('log')
plt.ylim(1,10)
plt.ylabel('Radio petrosiano 50$\%$')
plt.xlabel('Magnitud absoluta en r ($M_r$)')
plt.title('Diagrama de Magnitud Absoluta vs Radio ($50\%$)')
plt.plot()
plt.savefig('/mnt/sda2/extragalactica/2_practico/imagenes/mabs_vs_r50_bulgeydisco.png', dpi=300,bbox_inches='tight')
