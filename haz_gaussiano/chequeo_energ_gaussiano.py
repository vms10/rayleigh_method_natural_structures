
"""
Este script chequea el balance energetico para el caso de haz gaussiano.
Sirve tanto si el perfil es una funcion analitica o si viene de una imagen.
"""
#%% importo paquetes
import numpy as np
import cmath
from scipy.integrate import quad
import matplotlib.pyplot as plt
import time
# Lo que sigue solo descomentarlo si es necesario volver a cargar los datos de las soluciones y de los valores de los parametros
#%% hay que cargar la solucion para R
#npzfile_R=np.load(r"C:\Users... ruta del archivo.npz")
#x_R=npzfile_R["solucion"]
#alfas_R=npzfile_R["alfas"]
#%%  hay que cargar la solucion para T
#npzfile_T=np.load(r"C:\Users... ruta del archivo.npz")
#x_T=npzfile_T["solucion"]
#alfas=npzfile_T["alfas"]
#%%
#eje_y_2=[]
#for i in x_R:
#    eje_y_2.append(np.absolute(i)**2)
#%%
#eje_y=[]
#for i in x_T:
#    eje_y.append(np.absolute(i)**2)
#
#%% parametros constitutivos, modo de polarizacion y angulo de incidencia
## Necesita los siguientes valores
#epsilon_1=1
#mu_1=1
#epsilon_2=3
#mu_2=1
##tita=0
#tita=20*np.pi/180 # el valor hay que ponerlo en grados y lo convierte a radianes
#sigma_1=1 # =mu es TE, =epsilon es TM
#sigma_2=3
#
#
#long=1 #longitud de onda 
#
#k_0=2*np.pi/long
#alfa_0 = np.sin(tita)*k_0*((epsilon_1*mu_1)**(0.5))  



#%% CHEQUEO ENERGETICO 
from itertools import repeat
beta_alfa = lambda alpha,epsilon,mu: cmath.sqrt(epsilon*mu*(k_0**2)-alpha**2) # beta_alfa funciona para el medio 1 y 2 cambiando los epsilon y mu

betas_T=np.array(list(map(beta_alfa,alfas,repeat(epsilon_2),repeat(mu_2))))
betas_R=np.array(list(map(beta_alfa,alfas_R,repeat(epsilon_1),repeat(mu_1))))

#%%

T_abs_cuad=np.array(eje_y)
R_abs_cuad=np.array(eje_y_2)


#%%
puntos_T=betas_T*T_abs_cuad
int_T=np.trapz(puntos_T,x=alfas)
int_T=sigma_1*int_T.real/(beta_alfa(alfa_0,epsilon_1,mu_1)*sigma_2)
puntos_R=betas_R*R_abs_cuad
int_R=np.trapz(puntos_R,x=alfas_R)
int_R=int_R.real/(beta_alfa(alfa_0,epsilon_1,mu_1))
lado_der=int_R+int_T
#%%
#doble_v=4*a  #necesita el valor de a (tama;o rugosidad)
#desv=2/doble_v
#def A_alfa(alfa):
#    return np.exp(-((alfa-alfa_0)/desv)**2)/(2*np.pi*desv*np.pi**0.5)

#%%
A_alfass=[]
for i in alfas_R:
    A_alfass.append(np.absolute(A_alfa(i))**2)
    
A_alfass=np.array(A_alfass)
puntos_aa=betas_R*A_alfass
int_aa=np.trapz(puntos_aa,x=alfas_R)
int_aa=int_aa.real/beta_alfa(alfa_0,epsilon_1,mu_1)

#%%  error porcentual
e=100*np.absolute(lado_der/int_aa-1)

#%%
