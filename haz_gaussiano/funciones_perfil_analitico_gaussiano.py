"""
Este script tiene todas las funciones necesarias para correr la implementacion
tanto de la ecuacion para los campos reflejados como la de los campos transmitidos.
Lo mejor es ir corriendo celda por celda. 
Aca el perfil es introducido mediante una funcion analitica y la iluminacion es haz gaussiano.
"""
#%%
import numpy as np
import cmath
from scipy.integrate import quad
import matplotlib.pyplot as plt
import time
from itertools import repeat
from scipy import interpolate
#%% #parametros constitutivos, modo de polarizacion y angulo de incidencia
epsilon_1=1
mu_1=1
epsilon_2=3
mu_2=1
#tita=0
tita=30*np.pi/180 # el valor hay que ponerlo en grados y lo convierte a radianes
sigma_1=1 # =mu es TE, =epsilon es TM
sigma_2=1
#%%  perfil a partir de una funcion analitica
 #Tengo que dar las relaciones entre h,d,lambda
long=1 #longitud de onda 
d=2*long # "periodo" de la sinusoidal
h=0.02*long #altura
a=d  # tama;o total. La igualdad vale solo en el caso particular considerado

k_0=2*np.pi/long
alfa_0 = np.sin(tita)*k_0*((epsilon_1*mu_1)**(0.5))

 #%%   forma del perfil 
def g(x):
    return (1+np.cos(2*np.pi*x/d))*h/2

#%% # plot del perfil esquematico
equis=np.linspace(-a/2,a/2,500)
fig, ax = plt.subplots(figsize=(2*6.4, 2*4.8))

plt.xlabel(r'$x$', fontsize=18)

plt.ylabel(r'$g(x)$', fontsize=18)  
plt.plot(equis,g(equis))


#%%

def integrand_real(x, u, v):
    return (np.exp(-1j*(u*x+v*g(x)))).real
    
def integrand_imag(x, u, v):
    return (np.exp(-1j*(u*x+v*g(x)))).imag

def D(u,v):
    I_R = quad(integrand_real, -a/2, a/2, args=(u,v),limit=100)
    I_I = quad(integrand_imag, -a/2, a/2, args=(u,v),limit=100)
    return (I_R[0]+1j*I_I[0])/(2*np.pi)
#%%
doble_v=4*a
desv=2/doble_v
def A_alfa(alfa):
    return np.exp(-((alfa-alfa_0)/desv)**2)/(desv*np.pi**0.5)

#%%
beta_alfa = lambda alpha,epsilon,mu: cmath.sqrt(epsilon*mu*(k_0**2)-alpha**2) # beta_alfa funciona para el medio 1 y 2 cambiando los epsilon y mu

M_prima_alfa = lambda prima, alfa: ((1-sigma_2/sigma_1)*(prima*alfa+beta_alfa(alfa, epsilon_2, mu_2)*beta_alfa(prima,  epsilon_1, mu_1))+(k_0**2)*(sigma_2*epsilon_1*mu_1/sigma_1-epsilon_2*mu_2))/(beta_alfa(alfa,  epsilon_2, mu_2)-beta_alfa(prima, epsilon_1, mu_1))

def T_m(prima,alfa):
    return M_prima_alfa(prima,alfa)*(D(prima-alfa,beta_alfa(alfa,epsilon_2,mu_2)-beta_alfa(prima,epsilon_1,mu_1))-np.sin(a*(prima-alfa)/2)/((prima-alfa)*np.pi))   

v_T1 = np.vectorize(T_m)

def b_T(prima):
    return -2*sigma_2*beta_alfa(prima,epsilon_1,mu_1)*A_alfa(prima)/sigma_1

N_alfa_prima = lambda  alfa,prima: ((1-sigma_2/sigma_1)*(alfa*prima-beta_alfa(prima,epsilon_2, mu_2)*beta_alfa(alfa,  epsilon_1, mu_1))+(k_0**2)*(sigma_2*epsilon_1*mu_1/sigma_1-epsilon_2*mu_2))/(beta_alfa(prima, epsilon_2, mu_2)+beta_alfa(alfa, epsilon_1, mu_1))
#%%

def cacho_1_real(y, u):
    return (A_alfa(y)*N_alfa_prima(y,u)*np.sin(a*(u-y)/2)/(u-y)).real

def cacho_1_imag(y, u):
    return (A_alfa(y)*N_alfa_prima(y,u)*np.sin(a*(u-y)/2)/(u-y)).imag

    
#%%
def cacho_1(u):
    I_R = quad(cacho_1_real, -np.inf, np.inf , args=(u),limit=2000)
    I_I = quad(cacho_1_imag, -np.inf, np.inf, args=(u),limit=2000)
    return (I_R[0]+1j*I_I[0])/np.pi


#%% CHEQUEO CACHO1. Esto es para verificar que el integrador de la celda de arriba este funcionando bien.
# Esta comentado para no confundir.
#def FUNC(y,u):
#    return A_alfa(y)*N_alfa_prima(y,u)*np.sin(a*(u-y)/2)/(u-y)
#v_FUNC = np.vectorize(FUNC)
#equises=np.arange(-desv*8,desv*8,1/(10*a))
#yses=v_FUNC(equises,34)
#from scipy import integrate
#inted=integrate.simps(yses,x=equises)/np.pi


#%%
def cacho_2_real(x,y, u):
    return (A_alfa(y)*N_alfa_prima(y,u)*np.exp(-1j*((u-y)*x+(beta_alfa(u,epsilon_2,mu_2)+beta_alfa(y,epsilon_1,mu_1))*g(x)))).real

def cacho_2_imag(x,y, u):
    return (A_alfa(y)*N_alfa_prima(y,u)*np.exp(-1j*((u-y)*x+(beta_alfa(u,epsilon_2,mu_2)+beta_alfa(y,epsilon_1,mu_1))*g(x)))).imag

from scipy import integrate

def cacho_2(u):
    I_R = integrate.dblquad(cacho_2_real, -np.inf, np.inf,  -a/2,  a/2,args=[u])
    I_I = integrate.dblquad(cacho_2_imag, -np.inf, np.inf,  -a/2,  a/2,args=[u])
    return -(I_R[0]+1j*I_I[0])/(2*np.pi)

#%% #chequeo CACHO2 
#def FUNC(y,u):
#    return A_alfa(y)*N_alfa_prima(y,u)*D((u-y),(beta_alfa(u,epsilon_2,mu_2)-beta_alfa(y,epsilon_1,mu_1)))
#v_FUNC = np.vectorize(FUNC)
#equises=np.arange(-desv*8,desv*8,1/(10*a))
#yses=v_FUNC(equises,2)
#from scipy import integrate
#inted=-integrate.simps(yses,x=equises)



#%%
def R_m(prima,alfa):
    return M_prima_alfa(alfa,prima)*(D(prima-alfa,beta_alfa(prima,epsilon_2,mu_2)-beta_alfa(alfa,epsilon_1,mu_1))-np.sin(a*(prima-alfa)/2)/((prima-alfa)*np.pi))   

v_R1 = np.vectorize(R_m)

def b_R(prima):
    return A_alfa(prima)*(beta_alfa(prima,epsilon_2,mu_2)-sigma_2*beta_alfa(prima,epsilon_1,mu_1)/sigma_1)+cacho_1(prima)+cacho_2(prima)



#%%



