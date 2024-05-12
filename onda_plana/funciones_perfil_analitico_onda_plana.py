"""
Este script tiene todas las funciones necesarias para correr la implementacion
tanto de la ecuacion para los campos reflejados como la de los campos transmitidos.
Lo mejor es ir corriendo celda por celda. 
Aca el perfil es introducido mediante una funcion analitica y la iluminacion es onda plana.
"""
#%% importo paquetes
import numpy as np
import cmath
from scipy.integrate import quad
import matplotlib.pyplot as plt
import time
#%% parametros constitutivos, modo de polarizacion y angulo de incidencia
epsilon_1=1
mu_1=1
epsilon_2=3
mu_2=1
#tita=0
tita=20*np.pi/180 # el valor hay que ponerlo en grados y lo convierte a radianes
sigma_1=1 # =mu es TE, =epsilon es TM
sigma_2=1

#%%  perfil a partir de una funcion analitica
 #Tengo que dar las relaciones entre h,d,lambda
long=1 #longitud de onda. El parametro que elijo su valor numerico tiene que ser 1. Si no le ponen 1, hay que normalizar luego en las curvas.
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
    I_R = quad(integrand_real, -a/2, a/2, args=(u,v),limit=200) # limit es el numero maximo de ciclos que realiza para que la integral converja. Se le puede cambiar la precision con la que quiero que converja la integral. Por default es del orden 10^-9
    I_I = quad(integrand_imag, -a/2, a/2, args=(u,v),limit=200)
    return (I_R[0]+1j*I_I[0])/(2*np.pi)


#%%
beta_alfa = lambda alpha,epsilon,mu: cmath.sqrt(epsilon*mu*(k_0**2)-alpha**2) # beta_alfa funciona para el medio 1 y 2 cambiando los epsilon y mu
v_ba = np.vectorize(beta_alfa)
#%%

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



