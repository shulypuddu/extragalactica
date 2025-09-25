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

galaxias_brillante = df0[df0['abs_model_mag_r'] > -15]
#esta galaxia no se xq quedo tan FEA, asique la filtro
df = df0[df0['abs_model_mag_r'] < -15]


print(df.columns)
z= df['z']
M_r= df['abs_model_mag_r']
#%%
plt.figure(figsize=(10,6))
plt.scatter(z,M_r, s=1, color='teal', marker='.')
plt.xlabel('Redshift ($z$)')
plt.ylabel('Magnitud absoluta ($M_r$)')
plt.title('Diagrama de Redshift vs Magnitud Absoluta')
plt.gca().invert_yaxis()  
plt.plot()
plt.show()  

#%% ----------histogramas de indices de color en un mismo figure
fig, axs = plt.subplots(2, 1, figsize=(14, 18)) # 2 filas, 1 columna

# Primer gráfico: scatter de M_r vs z
axs[0].hist(df['color_g_r'], bins='auto', alpha=0.55,edgecolor='teal',facecolor='teal',density=True)
axs[0].set_xlabel('Índice de Color (g-r)')
axs[0].set_xlim(-0.25, 1.25)
axs[0].set_ylabel('Número de Galaxias')
axs[0].set_title('Histograma color_g_r')

# Segundo gráfico: histograma de color_u_r
axs[1].hist(df['color_u_r'], bins='auto', alpha=0.55,edgecolor='hotpink', facecolor='hotpink',density=True)
axs[1].set_xlabel('Índice de Color (u-r)')
axs[1].set_xlim(0, 3.65)
axs[1].set_ylabel('Número de Galaxias')
axs[1].set_title('Histograma color_u_r')

plt.tight_layout()
plt.show()

#Ahora voy a ajustar una doble gaussiana a cada figura
#%%
# Supón que tus datos están en la variable color_u_r
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
azul= gauss(x1_fit,params[0],params[1],params[2])
rojo= gauss(x2_fit,params[3],params[4],params[5])

# Graficar
plt.hist(color_u_r, bins=1000, density=True, alpha=0.5, edgecolor='gray', facecolor='lightgray', label='Datos')
plt.plot(x_fit, doble_gaussiana(x_fit, *params), label='Ajuste Doble Gaussiana', color='m')
plt.plot(x1_fit, azul,  label='Doble Gaussiana',color='skyblue') 
plt.plot(x2_fit, rojo, color='coral', label='Doble Gaussiana') 
plt.xlim(0, 4)
plt.xlabel('color_u_r') 
plt.ylabel('Densidad')
plt.show()

def diferencia(x):
    return gauss(x,params[0],params[1],params[2]) - gauss(x,params[3],params[4],params[5])
x_inicial= (params[4]+params[1])/2
x_cruce = opt.fsolve(diferencia, x_inicial)[0]

print(x_cruce)

#%%

fig, axs = plt.subplots(2, 1, figsize=(14, 18)) # 2 filas, 1 columna
# Primer gráfico: scatter de M_r vs z
axs[0].hist(df['color_g_r'], bins='auto',density=True, alpha=0.55,edgecolor='teal',facecolor='teal')
axs[0].set_xlabel('Índice de Color (g-r)')
axs[0].set_xlim(-0.6, 2)
axs[0].set_ylabel('Número de Galaxias')
axs[0].set_title('Histograma color_g_r')
# Segundo gráfico: histograma de color_u_r
axs[1].hist(df['color_u_r'], bins='auto',density=True, alpha=0.55,edgecolor='hotpink', facecolor='hotpink')
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
plt.xlim(1.15,4.0)
plt.title('Distribución del parámetro de concentración (C)')

#fracDeV habla de cuanto bulge hay en la galaxia.

plt.figure(figsize=(10,6))
plt.hist(df['fracDeV_r'], bins='auto',edgecolor='teal',facecolor='none', density=True)
plt.xlim(-0.0001,1.0001)
plt.title('Distribución del parámetro fracDeV en el filtro r')

plt.figure(figsize=(10,6))
plt.scatter(df['c_par'],df['fracDeV_r'], marker='.',alpha=0.5)
plt.xlabel('Parámetro de concentración (C)')
plt.ylabel('Frac de V')
plt.xlim(0.75,4)

#%% Creo que este es el que pide en el inciso donde habla de correlación
#Todo esto hay que hacerlo ya separando segun si son de la nube azul o de la secuencia roja 

#plt.title()

plt.figure(figsize=(10,6))
plt.scatter(df['color_u_r'],df['c_par'], marker='.',alpha=0.5)
plt.ylabel('(C)')
plt.xlabel('Indice de color u-r')
plt.xlim(0.75,4)



# DIAGRAMA COLOR MAGNITUD
plt.figure(figsize=(10,6))
plt.scatter(df['abs_petro_mag_r'], df['color_u_r'], marker='.',color='hotpink')
plt.xlim(-1,13)
plt.ylabel('Índice de color (u-r)')
plt.xlabel('$M_r$')
plt.title('Diagrama color - magnitud')
plt.gca().invert_xaxis()  # Invertir eje X (magnitud)
plt.plot()
plt.show()





#%% ------ Ahora separo por galaxias rojas y azules
rojo = df[df['color_u_r'] > x_cruce]
azul = df[df['color_u_r'] < x_cruce]

plt.figure(figsize=(10,6))
plt.scatter(rojo['color_u_r'],rojo['abs_petro_mag_r'],marker='.',color='coral',alpha=0.5)
plt.scatter(azul['color_u_r'],azul['abs_petro_mag_r'],marker='.',color='skyblue',alpha=0.5)
plt.xlim(-1,13)
plt.xlabel('Índice de color (u-r)')
plt.ylabel('$M_r$')
plt.title('Diagrama color - magnitud')
plt.plot()


#%%

plt.figure(figsize=(10,6))
plt.scatter(df['abs_model_mag_r'], df['petroR50_r'], s=1, color='teal', marker='.')
plt.ylabel('Radio que contiene el $50\%$ de la luz en r (petroR50_r)')
plt.xlabel('Magnitud absoluta en r ($M_r$)')
plt.title('Diagrama de Magnitud Absoluta vs Radio ($50\%$)')
plt.plot()

plt.figure(figsize=(10,6))
plt.scatter(rojo['abs_model_mag_r'], rojo['petroR50_r'], s=1, color='coral', marker='.',alpha=0.5)
plt.scatter(azul['abs_model_mag_r'], azul['petroR50_r'], s=1, color='skyblue', marker='.',alpha=0.5)
plt.ylabel('Radio que contiene el $50\%$ de la luz en r (petroR50_r)')
plt.xlabel('Magnitud absoluta en r ($M_r$)')
plt.title('Diagrama de Magnitud Absoluta vs Radio ($50\%$)')
plt.plot()

