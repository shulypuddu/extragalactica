#GRAFICOS PRACTICO 2 EXTRAGALAGTICA.

#%%
import pandas as pd
import matplotlib.pyplot as plt
b='royalblue'
rojo='tomato'
y='gold'
tb=pd.read_csv('pr2_Martina_Tetzlaff.csv')
z=tb['z']
ra=tb['ra']
dec=tb['dec']

plt.scatter(ra,dec,color=b,s=1)
plt.xlabel('Ascencion recta [°]')
plt.ylabel('declinación [°]')
plt.title('Distribución espacial de galaxias')
plt.show()


plt.hist(z,50,color=rojo)
plt.xlabel('redshift')
plt.title('Histograma de redshift')
plt.show()

r=tb['petroMag_r']

plt.scatter(z,r-tb['extinction_r'],color=b,s=1,label='Muestra reducida')
plt.xlabel('redshift')
plt.ylabel('r')

pt_r_c=[]
z_r=[]
for i in range(len(r)):
    if r[i]<14.5 or r[i]>17.77:
        pt_r_c.append(r[i]-tb['extinction_r'][i])
        z_r.append(z[i])
plt.scatter(z_r,pt_r_c,color=rojo,s=1,label='Galaxias descartadas')
plt.xlabel('redshift')
plt.ylabel('r')
plt.legend(loc='center right')
plt.title('Magnitud r aparente vs redshift')

#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

pr2_tb=pd.read_csv('tabla_pr2.csv')
#ID,ra,dec,petroM_u,petroM_g,petroM_r,petroM_i,petroM_z,modelM_u,modelM_g,modelM_r,modelM_i,modelM_z, u-r,g-r,C,fracDeV_r, mu_sup,velDisp,

com_vol=(pr2_tb['petroM_r']>-20.3) & (pr2_tb['petroM_r']<-19)
tb_m=pr2_tb

plt.scatter(tb_m['z'],tb_m['petroM_r'],color=b,s=3,label='Muestra completa por flujo')
plt.xlabel('z')
plt.ylabel('$M_r$')
plt.gca().invert_yaxis()
plt.yscale('linear')
plt.legend()
plt.title('Magnitud r vs redshift')
plt.show()


#-19.1,-20.3

from scipy.optimize import curve_fit
from scipy.optimize import fsolve

bin=50
plt.hist(tb_m['u-r'],range=[0.0,3.5],bins=bin,density=True,color='k',alpha=0.4,label='Muestra')

def gausss(x, A1, mu1, sigma1, A2, mu2, sigma2):
    gauss1 = A1 * np.exp(-(x - mu1)**2 / (2*sigma1**2))
    gauss2 = A2 * np.exp(-(x - mu2)**2 / (2*sigma2**2))
    return gauss1 + gauss2


x=np.linspace(0,4,bin)
counts,edge = np.histogram(tb_m['u-r'], range=[0.0, 4], bins=bin, density=True)
p0 = [2500, 1.5, 0.2, 2300, 2.5, 0.1]   
popt, pcov = curve_fit(gausss, x, counts, p0=p0)

def gauss(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2*sigma**2))


plt.plot(x,gauss(x,popt[0],popt[1],popt[2]),color=b,label='Galaxias azules')
plt.plot(x,gauss(x,popt[3],popt[4],popt[5]),color=rojo,label='Galaxias rojas')
plt.xlabel('u-r')
plt.title('Histograma color u-r')
plt.legend()
plt.show()

def h(y):
    h=gauss(y,popt[0],popt[1],popt[2])-gauss(y,popt[3],popt[4],popt[5])
    return h
x0 = 2.3  
lim_color = fsolve(h, x0)
print(lim_color)

#------------------------------------------------------

x=np.linspace(0,1.3,bin)
counts,edge = np.histogram(tb_m['g-r'], range=[0.0, 1.3], bins=bin)
p0 = [10000, 0.4, 0.1, 10000, 0.8, 0.5]   
popt, pcov = curve_fit(gausss, x, counts, p0=p0)

lim_color = float(lim_color)
mask = tb_m['u-r'] <= lim_color   
g_azul = tb_m[mask]
g_rojo = tb_m[~mask]

plt.hist(tb_m['g-r'],range=[0.0,1.3],bins=bin,color='k',alpha=0.4,label='Muestra')
plt.plot(x,gausss(x,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]),color='goldenrod',linewidth=2.2,label='Ajuste bimodal')
plt.hist(g_azul['g-r'],range=[0.0,1.3],bins=bin,color=b,alpha=0.3,label='Galaxias azules')
#plt.plot(x,gauss(x,popt[3],popt[4],popt[5]),color=rojo,label='Galaxias rojas')
plt.hist(g_rojo['g-r'],range=[0.0,1.3],bins=bin,color=rojo,alpha=0.3,label='Galaxias rojas')

plt.xlabel('g-r')
plt.title('Histograma color g-r')
plt.legend()
plt.show()

#%%
import seaborn as sns

sns.jointplot(data=tb_m, x='C', y='fracDeV_r', kind='scatter', color=b,s=3, marginal_kws=dict(bins='auto', fill=True))
plt.axvline(2.5, color=rojo, linestyle='solid', label='Valor crítico', linewidth=3)
plt.xlim(1.3, 4)
plt.ylim(0, 1)
plt.show()
#%%
#Parametro de concentración vs color
import matplotlib.patches as patches

fig, ax = plt.subplots()
plt.scatter(pr2_tb['C'],pr2_tb['u-r'],color='k',alpha=0.4,s=0.5)
plt.ylim((0.5,3.5))
plt.xlim(1.3,4)
plt.title('Diagrama parametro de concentración vs. color')
plt.xlabel('C')
plt.ylabel('u-r')
plt.vlines(2.5,0.5,3.5,'k','dashed')
plt.hlines(lim_color,1.3,4,'k','dashed')
rect = patches.Rectangle((1.3, 0.5), 1.2, lim_color-0.5, linewidth=1, edgecolor='k', facecolor=b, alpha=0.2, label='Azules y tardías')
ax.add_patch(rect)
rect = patches.Rectangle((2.5, lim_color), 4, 3.5, linewidth=1, edgecolor='k', facecolor=rojo, alpha=0.2, label='Rojas y tempranas')
ax.add_patch(rect)
plt.legend()
plt.show



#%% 
#JULIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII ACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
m_25=np.percentile(tb_m['petroM_r'],25)
mask_25=tb_m['petroM_r']<=m_25
m_50=np.percentile(tb_m['petroM_r'],50)
mask_50=(tb_m['petroM_r']<=m_50) & (tb_m['petroM_r']>m_25)
m_75=np.percentile(tb_m['petroM_r'],75)
mask_75=(tb_m['petroM_r']<=m_75) & (tb_m['petroM_r']>m_50)
mask_100=tb_m['petroM_r']>m_75

from scipy.optimize import curve_fit
from scipy.optimize import fsolve

bin=50
plt.hist(tb_m['u-r'][mask_25],range=[0.0,3.5],bins=bin,density=True,color='k',alpha=0.4,label=f'Muestra $M_r$ = [-22, {m_25:.2f}]')
def gausss(x, A1, mu1, sigma1, A2, mu2, sigma2):
    gauss1 = A1 * np.exp(-(x - mu1)**2 / (2*sigma1**2))
    gauss2 = A2 * np.exp(-(x - mu2)**2 / (2*sigma2**2))
    return gauss1 + gauss2


x=np.linspace(0,4,bin)
counts,edge = np.histogram(tb_m['u-r'][mask_25], range=[0.0, 4], bins=bin, density=True)
p0 = [0.8, 1.9, 0.2, 1.6, 2.6, 0.1]   
popt, pcov = curve_fit(gausss, x, counts, p0=p0)

def gauss(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2*sigma**2))

mu_25a=popt[1]
sig_25a=popt[2]
mu_25r=popt[4]
sig_25r=popt[5]
plt.plot(x,gauss(x,popt[0],popt[1],popt[2]),color=b,label='Galaxias azules')
plt.plot(x,gauss(x,popt[3],popt[4],popt[5]),color=rojo,label='Galaxias rojas')
plt.xlabel('u-r')
plt.title('Histograma color u-r primer cuartil')
plt.legend()
plt.show()
#----------------------------------------------------------------------------
bin=50
plt.hist(tb_m['u-r'][mask_50],range=[0.0,3.5],bins=bin,density=True,color='k',alpha=0.4,label=f'Muestra $M_r$ = [{m_25:.2f}, {m_50:.2f}]')
def gausss(x, A1, mu1, sigma1, A2, mu2, sigma2):
    gauss1 = A1 * np.exp(-(x - mu1)**2 / (2*sigma1**2))
    gauss2 = A2 * np.exp(-(x - mu2)**2 / (2*sigma2**2))
    return gauss1 + gauss2


x=np.linspace(0,4,bin)
counts,edge = np.histogram(tb_m['u-r'][mask_50], range=[0.0, 4], bins=bin, density=True)
p0 = [0.8, 1.9, 0.2, 1.6, 2.6, 0.1]   
popt, pcov = curve_fit(gausss, x, counts, p0=p0)

def gauss(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2*sigma**2))

mu_50a=popt[1]
sig_50a=popt[2]
mu_50r=popt[4]
sig_50r=popt[5]
plt.plot(x,gauss(x,popt[0],popt[1],popt[2]),color=b,label='Galaxias azules')
plt.plot(x,gauss(x,popt[3],popt[4],popt[5]),color=rojo,label='Galaxias rojas')
plt.xlabel('u-r')
plt.title('Histograma color u-r segundo cuartil')
plt.legend()
plt.show()

#----------------------------------------------------------------------------
bin=50
plt.hist(tb_m['u-r'][mask_75],range=[0.0,3.5],bins=bin,density=True,color='k',alpha=0.4,label=f'Muestra $M_r$ = [{m_50:.2f}, {m_75:.2f}]')
def gausss(x, A1, mu1, sigma1, A2, mu2, sigma2):
    gauss1 = A1 * np.exp(-(x - mu1)**2 / (2*sigma1**2))
    gauss2 = A2 * np.exp(-(x - mu2)**2 / (2*sigma2**2))
    return gauss1 + gauss2


x=np.linspace(0,4,bin)
counts,edge = np.histogram(tb_m['u-r'][mask_75], range=[0.0, 4], bins=bin, density=True)
p0 = [0.8, 1.9, 0.2, 1.6, 2.6, 0.1]   
popt, pcov = curve_fit(gausss, x, counts, p0=p0)

def gauss(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2*sigma**2))

mu_75a=popt[1]
sig_75a=popt[2]
mu_75r=popt[4]
sig_75r=popt[5]
plt.plot(x,gauss(x,popt[0],popt[1],popt[2]),color=b,label='Galaxias azules')
plt.plot(x,gauss(x,popt[3],popt[4],popt[5]),color=rojo,label='Galaxias rojas')
plt.xlabel('u-r')
plt.title('Histograma color u-r tercer cuartil')
plt.legend()
plt.show()

bin=50
plt.hist(tb_m['u-r'][mask_100],range=[0.0,3.5],bins=bin,density=True,color='k',alpha=0.4,label=f'Muestra $M_r$ = [{m_75:.2f}, -17]')
def gausss(x, A1, mu1, sigma1, A2, mu2, sigma2):
    gauss1 = A1 * np.exp(-(x - mu1)**2 / (2*sigma1**2))
    gauss2 = A2 * np.exp(-(x - mu2)**2 / (2*sigma2**2))
    return gauss1 + gauss2


x=np.linspace(0,4,bin)
counts,edge = np.histogram(tb_m['u-r'][mask_100], range=[0.0, 4], bins=bin, density=True)
p0 = [0.8, 1.9, 0.2, 1.6, 2.6, 0.1]   
popt, pcov = curve_fit(gausss, x, counts, p0=p0)

def gauss(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2*sigma**2))

mu_100a=popt[1]
sig_100a=popt[2]
mu_100r=popt[4]
sig_100r=popt[5]
plt.plot(x,gauss(x,popt[0],popt[1],popt[2]),color=b,label='Galaxias azules')
plt.plot(x,gauss(x,popt[3],popt[4],popt[5]),color=rojo,label='Galaxias rojas')
plt.xlabel('u-r')
plt.title('Histograma color u-r cuarto cuartil')
plt.legend()
plt.show()

#%%
#%%
#Magnitud vs color separado por azules y rojas



# grafico
mas_c=tb_m['C']>2.5
import matplotlib.pyplot as plt
plt.figure(figsize=(8,6))
plt.scatter(tb_m['petroM_r'], tb_m['u-r'], color=b, label='Tardías',s=1)
plt.scatter(tb_m['petroM_r'][mas_c], tb_m['u-r'][mas_c], color=rojo, label='Tempranas',s=1)
plt.vlines((m_25-22)/2,mu_25a-sig_25a,mu_25a+sig_25a,'b','solid',alpha=0.5,linewidth=8,label='Dispersión de las galaxias azules')
plt.vlines((m_50+m_25)/2,mu_50a-sig_50a,mu_50a+sig_50a,'b','solid',alpha=0.5,linewidth=8)
plt.vlines((m_75+m_50)/2,mu_75a-sig_75a,mu_75a+sig_75a,'b','solid',alpha=0.5,linewidth=8)
plt.vlines((-17+m_75)/2,mu_100a-sig_100a,mu_50a+sig_50a,'b','solid',alpha=0.5,linewidth=8)

plt.vlines((m_25-22)/2,mu_25r-sig_25r,mu_25r+sig_25r,'r','solid',alpha=0.5,linewidth=8,label='Dispersión de las galaxias rojas')
plt.vlines((m_50+m_25)/2,mu_50r-sig_50r,mu_50r+sig_50r,'r','solid',alpha=0.5,linewidth=8)
plt.vlines((m_75+m_50)/2,mu_75r-sig_75r,mu_75r+sig_75r,'r','solid',alpha=0.5,linewidth=8)
plt.vlines((-17+m_75)/2,mu_100r-sig_100r,mu_50r+sig_50r,'r','solid',alpha=0.5,linewidth=8)

plt.xlabel('$M_r$')
plt.ylabel('u-r')
plt.gca().invert_xaxis()
plt.ylim((0.5,3))
plt.hlines(lim_color,-22,-17,'k','dashed',label='Valor crítico de u-r')
plt.title('Diagrama color-magnitud')
plt.legend()
plt.show()

def h(x,a,b):
    h=a*x+b
    return h
p0 = [-1,3]   
popt, pcov = curve_fit(h, g_azul['petroM_r'], g_azul['u-r'], p0=p0)

x=np.linspace(-22,-17,50)

#dividir en cuatriles para M asi ajustamos gaussianas a u-r por secciones.


#%%
#Relacion de cormendy log y mu sup, proyeccion fotometrica (porque no tiene un cuenta la dispersion de velocidades-) del plano fundamental donde viven las calaxias espirales.
lim_color = float(lim_color)   # forzar a escalar
pr2_tb=pd.read_csv('opt_tablaprac2.csv')
mask = pr2_tb['u-r'] <= lim_color   # mask es una Series booleana con la misma longitud que pr2_tb
grupo_azul = pr2_tb[mask]
grupo_rojo = pr2_tb[~mask]

coeffs = np.polyfit(np.log(grupo_azul['R50']), grupo_azul['petroM_r'], 1)  # 1 means linear
a_a, b_a = coeffs
coeffs = np.polyfit(np.log(grupo_rojo['R50']), grupo_rojo['petroM_r'], 1)  # 1 means linear
a_r, b_r = coeffs
x=np.linspace(0.6,9,2)
# grafico
import matplotlib.pyplot as plt
plt.figure(figsize=(8,6))
plt.scatter(grupo_azul['R50'], grupo_azul['petroM_r'], color=b, label='Galaxias azules', alpha=0.5,s=3)
plt.scatter(grupo_rojo['R50'], grupo_rojo['petroM_r'], color=rojo, label='Galaxias rojas', alpha=0.5,s=3)
plt.plot(x,a_a*np.log(x)+b_a,color='b',label='Ajuste lineal azul')
plt.plot(x,a_r*np.log(x)+b_r,color='r',label='Ajuste lineal rojo')
plt.xlabel('R50')
plt.xscale('log')
plt.xlim(0.55,10)
plt.gca().invert_yaxis()
plt.legend()
plt.ylabel('$M_r$')
plt.show()

#%%
#Relacion de cormendy log y mu sup, proyeccion fotometrica (porque no tiene un cuenta la dispersion de velocidades-) del plano fundamental donde viven las calaxias espirales.
lim_color = float(lim_color)   # forzar a escalar
pr2_tb=pd.read_csv('opt_tablaprac2.csv')
mask = pr2_tb['u-r'] <= lim_color   # mask es una Series booleana con la misma longitud que pr2_tb
grupo_azul = pr2_tb[mask]
grupo_rojo = pr2_tb[~mask]

coeffs = np.polyfit(np.log(grupo_azul['R50'])[grupo_azul['C']<2.5], grupo_azul['petroM_r'][grupo_azul['C']<2.5], 1)  # 1 means linear
a_a, b_a = coeffs
coeffs = np.polyfit(np.log(grupo_rojo['R50'])[grupo_rojo['C']>2.5], grupo_rojo['petroM_r'][grupo_rojo['C']>2.5], 1)  # 1 means linear
a_r, b_r = coeffs
x=np.linspace(0.6,9,2)

# grafico
import matplotlib.pyplot as plt
plt.figure(figsize=(8,6))
plt.scatter(grupo_azul['R50'][grupo_azul['C']<2.5], grupo_azul['petroM_r'][grupo_azul['C']<2.5], color=b, label='Galaxias azules y tardías', alpha=0.5,s=3)
plt.scatter(grupo_rojo['R50'][grupo_rojo['C']>2.5], grupo_rojo['petroM_r'][grupo_rojo['C']>2.5], color=rojo, label='Galaxias rojas y tempranas', alpha=0.5,s=3)
plt.xlabel('R50')
plt.xscale('log')
plt.plot(x,a_a*np.log(x)+b_a,color='b',label='Ajuste lineal azul')
plt.plot(x,a_r*np.log(x)+b_r,color='r',label='Ajuste lineal rojo')
plt.legend()
plt.gca().invert_yaxis()
plt.xlim(0.55,10)
plt.ylabel('$M_r$')
plt.show()


#%%DESCARTADOSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

#%%

lim_color = float(lim_color)   # forzar a escalar
pr2_tb=pd.read_csv('opt_tablaprac2.csv')
mask = pr2_tb['u-r'] <= lim_color   # mask es una Series booleana con la misma longitud que pr2_tb
grupo_azul = pr2_tb[mask]
grupo_rojo = pr2_tb[~mask]

# grafico
import matplotlib.pyplot as plt
plt.figure(figsize=(8,6))
plt.scatter(grupo_azul['R50'], grupo_azul['petroM_r'], color='b', label=f'u-r ≤ {lim_color}', alpha=0.6)
plt.scatter(grupo_rojo['R50'], grupo_rojo['petroM_r'], color='r', label=f'u-r > {lim_color}', alpha=0.6)
plt.xlabel('R50')
plt.xscale('log')
plt.ylabel('petroM_r')
plt.show()

#Hacer este grafico con las otras subdivisiones.


#%%
from scipy.optimize import curve_fit

bin=50
pr2_tb=pd.read_csv('opt_tablaprac2.csv')
plt.hist(pr2_tb['g-r'],range=[0.0,1.3],bins=bin,density=True,color='darkseagreen')

def gausss(x, A1, mu1, sigma1, A2, mu2, sigma2):
    gauss1 = A1 * np.exp(-(x - mu1)**2 / (2*sigma1**2))
    gauss2 = A2 * np.exp(-(x - mu2)**2 / (2*sigma2**2))
    return gauss1 + gauss2


x=np.linspace(0,1.3,bin)
counts,edge = np.histogram(pr2_tb['g-r'], range=[0.0, 1.3], bins=bin, density=True)
p0 = [0.7, 0.4, 0.1, 0.8, 0.8, 0.05]   
popt, pcov = curve_fit(gausss, x, counts, p0=p0)

def gauss(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2*sigma**2))

plt.plot(x,gauss(x,popt[0],popt[1],popt[2]),color='purple')
plt.plot(x,gauss(x,popt[3],popt[4],popt[5]),color='orange')
plt.xlabel('g-r')
plt.title('Histograma color g-r')


#%%
#Parametro de concentración vs color
plt.scatter(pr2_tb['C'],pr2_tb['u-r'],color='seagreen',marker='.',alpha=0.3)
plt.ylim((0,4))
plt.xlim((0,4))
plt.title('Diagrama parametro de concentración vs. color')
plt.show


#%%
#Magnitud vs color separado por azules y rojas

lim_color = float(lim_color)   # forzar a escalar
pr2_tb=pd.read_csv('opt_tablaprac2.csv')
mask = pr2_tb['u-r'] <= lim_color   
grupo_azul = pr2_tb[mask]
grupo_rojo = pr2_tb[~mask]

# grafico
import matplotlib.pyplot as plt
plt.figure(figsize=(8,6))
plt.scatter(grupo_azul['petroM_r'], grupo_azul['u-r'], color=b, label=f'u-r ≤ {lim_color}',s=2,alpha=0.7)
plt.scatter(grupo_rojo['petroM_r'], grupo_rojo['u-r'], color=rojo, label=f'u-r > {lim_color}',s=2,alpha=0.7)
plt.xlabel('$M_r$')
plt.ylabel('u-r')
plt.gca().invert_xaxis()
plt.ylim((0.5,3))
plt.title('Magnitud r vs color')
plt.show()

def h(x,a,b):
    h=a*x+b
    return h
p0 = [-1,3]   
popt, pcov = curve_fit(h, grupo_azul['petroM_r'], grupo_azul['u-r'], p0=p0)

x=np.linspace(-22,-17,50)

#dividir en cuatriles para M asi ajustamos gaussianas a u-r por secciones.

# %%
#Histograma de parametro de concentración
plt.hist(tb_m['C'],density=True,range=[1.5,4],bins='auto',color=b)
plt.title('Histograma parámetro de concentración')
plt.xlabel('')
plt.vlines(2.5,0,1.05,'darkblue','solid',label='Valor crítico')
plt.legend()
plt.show()

#valor estandar 2.5 para cortar entre tempranas y tardias
#Histpgrama de fracción
plt.hist(tb_m['fracDeV_r'],density=True,bins='auto',color=rojo)
plt.title('Histograma de la fracción DeVaculers')
plt.xlabel('fracDeV_r')
plt.show()

#Parametro de concentración vs francción
plt.scatter(tb_m['C'],tb_m['fracDeV_r'],color=rojo,marker='.',alpha=0.3)
plt.title('Parámetro de concentración vs. fracción DeVaculers')
plt.xlabel('C')
plt.ylabel('fracDeV_r')
plt.vlines(2.5,0,1,b,'solid',label='Valor crítico de c',linewidth=3)
plt.xlim(1.3,4)
plt.legend()
plt.show()
#Vacolours tardio, exp temprana

#%%
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(8, 8))
gs = gridspec.GridSpec(2, 2, width_ratios=[4, 1], height_ratios=[1, 4],
                       wspace=0.05, hspace=0.05)

# Main scatter plot
ax_main = plt.subplot(gs[1, 0])
ax_main.scatter(tb_m['C'], tb_m['fracDeV_r'], color=rojo, marker='.', alpha=0.3)
ax_main.vlines(2.5, 0, 1, b, 'solid', label='Valor crítico de c', linewidth=3)
ax_main.set_xlim(1.3, 4)
ax_main.set_xlabel('C')
ax_main.set_ylabel('fracDeV_r')
ax_main.legend()
ax_main.set_title('Parámetro de concentración vs. fracción DeVaculers')

# Top histogram (C)
ax_top = plt.subplot(gs[0, 0], sharex=ax_main)
ax_top.hist(tb_m['C'], density=True, range=[1.5, 4], bins='auto', color=b)
ax_top.vlines(2.5, 0, 1.05, 'darkblue', 'solid', label='Valor crítico')
ax_top.axis('off')  # Hide axes for cleaner look

# Right histogram (fracDeV_r)
ax_right = plt.subplot(gs[1, 1], sharey=ax_main)
ax_right.hist(tb_m['fracDeV_r'], density=True, bins='auto', color=rojo, orientation='horizontal')
ax_right.axis('off')  # Hide axes for cleaner look

plt.show()