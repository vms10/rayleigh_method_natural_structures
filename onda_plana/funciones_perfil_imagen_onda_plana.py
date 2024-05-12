"""
Este script tiene todas las funciones necesarias para correr la implementacion
tanto de la ecuacion para los campos reflejados como la de los campos transmitidos.
Lo mejor es ir corriendo celda por celda. 
Aca el perfil es introducido mediante coordenadas de una imagen y la iluminacion es onda plana.
"""
#%% importo paquetes
import numpy as np
import cmath
from scipy.integrate import quad
import matplotlib.pyplot as plt
import time
from itertools import repeat
from scipy import interpolate
#%% parametros constitutivos, modo de polarizacion y angulo de incidencia
epsilon_1=1.35**2
mu_1=1
epsilon_2=1.5**2
mu_2=1
tita=0
#tita=20*np.pi/180 # el valor hay que ponerlo en grados y lo convierte a radianes
sigma_1=1 # =mu es TE, =epsilon es TM
sigma_2=1
#%%  funcion para importar el txt con las coordenadas del perfil
def ldata(archive):
    f=open(archive)
    data=[]
    for line in f:
        e=line.strip()
        col=line.split()
        data.append(col)	
    return data

    #%%  importo el txt
    
data=ldata("euglenido_paper_listo.txt")
xfin=[]
for i in data:
    xfin.append(float(i[0]))
    
yfin=[]
for i in data:
    yfin.append(float(i[1]))
  #%%
# verificar que  1) xfin[0]=-xfin[-1]   # simetrico respecto de x=0
              #  2) yfin[0]=yfin[-1]=0  #en los bordes valga g(x=+/-a/2)=0
#%%  
g = interpolate.interp1d(xfin, yfin)
xneww = np.linspace(xfin[0], xfin[-1], 10**4) # 10**4 se puede cambiar
yneww = g(xneww)   # use interpolation function returned by `interp1d`


from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
fig, ax = plt.subplots(figsize=(3*6.50,3*1.38))

ax.plot(xfin,yfin,'ro',markersize=2)
ax.plot(xneww,yneww,'b-',markersize=2)
#ax.scatter(x_cuatro,y_cuatro)  
plt.xlabel(r'$x(\mu m )$', fontsize=24, labelpad=20)
ax.yaxis.set_major_locator(plt.MaxNLocator(3))
ax.tick_params(which='both', width=2,length=7, color='k')
plt.ylabel(r'$y(\mu m )$', fontsize=24, labelpad=20)
#plt.ylim(-0.005,0.02)
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(20) 
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(20) 
#ax.set_xticklabels(fontsize=12) 
#ax.annotate(r'$\nu_{1}=1.35 $',(-0.65,0.0125),fontsize=24)
#ax.annotate(r'$\nu_{2}=1.5 $',(-0.65,-0.003),fontsize=24)
#ax.xaxis.set_minor_locator(MultipleLocator(0.1))
#ax.yaxis.set_minor_locator(MultipleLocator(0.005))
plt.tight_layout() 

#%% 

a=xfin[-1]-xfin[0] # tama;o de la rugosidad
long=0.3 #longitud de onda en micrones
k_0=2*np.pi/long
alfa_0 = np.sin(tita)*k_0*((epsilon_1*mu_1)**(0.5))


#%% integrador 
def integrand(x, y, u, v):
    return np.exp(-1j*(u*x+v*y))
    
v_integrand = np.vectorize(integrand)
  
from scipy import integrate

def Dm(xnew,ynew,u,v):
    return integrate.simps(v_integrand(xnew,ynew,u,v), x=xnew)/(2*np.pi)

#%%
def D(u,v):
    N = 70 #70
    xnew = np.linspace(xfin[0], xfin[-1], N) # 10**4 se puede cambiar
    ynew = g(xnew)   # use interpolation function returned by `interp1d`
    
    old_integral = Dm(xnew,ynew,u,v)
    while True:
        N *= 2 #2
        xnew = np.linspace(xfin[0], xfin[-1], N) # 10**4 se puede cambiar
        ynew = g(xnew)   # use interpolation function returned by `interp1d`
        new_integral = Dm(xnew,ynew,u,v)
        if np.abs(old_integral - new_integral) < 1.49e-8:
            return (4*new_integral - old_integral)/3
        old_integral = new_integral

#%%
beta_alfa = lambda alpha,epsilon,mu: cmath.sqrt(epsilon*mu*(k_0**2)-alpha**2) # beta_alfa funciona para el medio 1 y 2 cambiando los epsilon y mu
v_ba = np.vectorize(beta_alfa)


M_prima_alfa = lambda prima, alfa: ((1-sigma_2/sigma_1)*(prima*alfa+beta_alfa(alfa, epsilon_2, mu_2)*beta_alfa(prima,  epsilon_1, mu_1))+(k_0**2)*(sigma_2*epsilon_1*mu_1/sigma_1-epsilon_2*mu_2))/(beta_alfa(alfa,  epsilon_2, mu_2)-beta_alfa(prima, epsilon_1, mu_1))

#%%  Funciones para la ecuacion de T
T = (2*sigma_2*beta_alfa(alfa_0,epsilon_1,mu_1)/sigma_1)/(sigma_2*beta_alfa(alfa_0,epsilon_1,mu_1)/sigma_1+beta_alfa(alfa_0,epsilon_2,mu_2)) # coef de fresnel

def T_m(alfa,prima):
    return M_prima_alfa(alfa,prima)*(D(alfa-prima,beta_alfa(prima,epsilon_2,mu_2)-beta_alfa(alfa,epsilon_1,mu_1))-np.sin(a*(alfa-prima)/2)/((alfa-prima)*np.pi))   

v_T1 = np.vectorize(T_m)

def b_T(alfa):
    return 2*T*M_prima_alfa(alfa,alfa_0)*(np.sin(a*(alfa-alfa_0)/2)/(alfa-alfa_0)-np.pi*D(alfa-alfa_0,beta_alfa(alfa_0,epsilon_2,mu_2)-beta_alfa(alfa,epsilon_1,mu_1)))


#%%   Funciones para la ecuacion de R
R = (sigma_2*beta_alfa(alfa_0,epsilon_1,mu_1)/sigma_1-beta_alfa(alfa_0,epsilon_2,mu_2))/(sigma_2*beta_alfa(alfa_0,epsilon_1,mu_1)/sigma_1+beta_alfa(alfa_0,epsilon_2,mu_2))  #coef de fresnel

N_alfa0_alfa = lambda  alfa: ((1-sigma_2/sigma_1)*(alfa_0*alfa-beta_alfa(alfa,epsilon_2, mu_2)*beta_alfa(alfa_0,  epsilon_1, mu_1))+(k_0**2)*(sigma_2*epsilon_1*mu_1/sigma_1-epsilon_2*mu_2))/(beta_alfa(alfa, epsilon_2, mu_2)+beta_alfa(alfa_0, epsilon_1, mu_1))


def R_m(alfa,prima):
    return M_prima_alfa(prima,alfa)*(D(alfa-prima,beta_alfa(alfa,epsilon_2,mu_2)-beta_alfa(prima,epsilon_1,mu_1))-np.sin(a*(alfa-prima)/2)/((alfa-prima)*np.pi))   

v_R1 = np.vectorize(R_m)

def b_R(prima):
    return (R*M_prima_alfa(alfa_0,prima)+N_alfa0_alfa(prima))*(2*np.sin(a*(prima-alfa_0)/2)/(prima-alfa_0))-2*np.pi*N_alfa0_alfa(prima)*D(prima-alfa_0,beta_alfa(prima,epsilon_2,mu_2)+beta_alfa(alfa_0,epsilon_1,mu_1))-2*np.pi*R*M_prima_alfa(alfa_0,prima)*D(prima-alfa_0,beta_alfa(prima,epsilon_2,mu_2)-beta_alfa(alfa_0,epsilon_1,mu_1))


#%%


