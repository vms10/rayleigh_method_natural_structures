
import numpy as np
import cmath
from scipy.integrate import quad
import matplotlib.pyplot as plt
import time
from itertools import repeat
from scipy import interpolate
#%%
def ldata(archive):
    f=open(archive)
    data=[]
    for line in f:
        e=line.strip()
        col=line.split()
        data.append(col)	
    return data
#%%  # txt con coordenadas del imagej. NO TIENE QUE TENER ENCABEZADOS
data=ldata("euglenidos_raw_paper2.txt")
xfin=[]
for i in data:
    xfin.append(float(i[0]))
    
yfin=[]
for i in data:
    yfin.append(float(i[1]))


 #%% # Plot para ver como se ven los puntos que exporte del ImageJ
fig, ax = plt.subplots(figsize=(2*6.4, 2*4.8))
ax.plot(xfin,yfin,"r", markersize=3)
plt.xlabel(r'$x(\mu m )$', fontsize=18)

plt.ylabel(r'$y(\mu m )$', fontsize=18)   
"""
Como se puede ver, ademas del perfil deseado hay lineas que no se quieren.
Para eliminarlas hay que explorar la imagen y delimitar los limites que
se quieren borrar. Estos bordes indeseados casi siempre aparecen y hay que eliminarlos. Los limites los identifique en el plot siguiente con
lineas horizontales y verticales:

"""

 #%% # Plot para ver bordes
fig, ax = plt.subplots(figsize=(2*6.4, 2*4.8))
ax.plot(xfin,yfin,"r" ,markersize=3)
plt.xlabel(r'$x(\mu m )$', fontsize=18)

plt.ylabel(r'$y(\mu m )$', fontsize=18)
ax.axvline( x=0.259)
ax.axvline( x=1.689)

ax.axhline( y=0.9)


#%% # Elimino bordes que no quiero en x
x_posta2=[]
y_posta2=[]
for i in range(len(xfin)):
    if (xfin[i]>0.259) & (xfin[i]<1.689):
        y_posta2.append(yfin[i])
        x_posta2.append(xfin[i])
    else:
        pass


#%% # Veo como quedaron los puntos
fig, ax = plt.subplots(figsize=(2*6.4, 2*4.8))
ax.plot(x_posta2,y_posta2,"ro", markersize=3)
plt.xlabel(r'$x(\mu m )$', fontsize=18)

plt.ylabel(r'$y(\mu m )$', fontsize=18)
ax.axhline( y=0.87)



#%% #elimino bordes en y
x_posta3=[]
y_posta3=[]
for i in range(len(x_posta2)):
    if (y_posta2[i]>0.87):
        x_posta3.append(x_posta2[i])
        y_posta3.append(y_posta2[i])
       

#%% #
fig, ax = plt.subplots(figsize=(2*6.4, 2*4.8))

ax.plot(x_posta3,y_posta3,"ro",markersize=2)
plt.xlabel(r'$x(\mu m )$', fontsize=18)

plt.ylabel(r'$y(\mu m )$', fontsize=18) 

ax.axvline( x=0.29)
ax.axhline( y=0.9)
    

plt.ylim(0,1.05)

  
#%% #elimino mas bordes  
x_posta4=[]
y_posta4=[]
for i in range(len(x_posta3)):
    if (x_posta3[i]<0.29) & (y_posta3[i]<0.9):
        pass
    else:
        y_posta4.append(y_posta3[i])
        x_posta4.append(x_posta3[i])

  #%% # Puntos finales sin ningun borde que no quiera 
fig, ax = plt.subplots(figsize=(2*6.4, 2*4.8))

ax.plot(x_posta4,y_posta4,"ro",markersize=2)
plt.xlabel(r'$x(\mu m )$', fontsize=18)

plt.ylabel(r'$y(\mu m )$', fontsize=18)    
plt.xlim(0.2,1.8)
plt.ylim(0,1.05)



#%% funcion para rotar y trasladar en las condiciones del codigo 
def rotar_trasladar (x,y): # x e y tienen que ser listas
    points=[]
    for i in range(len(x)):
        points.append([x[i],y[i]])
    points=sorted(points)    
 #El angulo que hay que rotar los puntos
    angulo=np.arctan((points[-1][1]-points[0][1])/(points[-1][0]-points[0][0]))
     # Matriz de rotacion 
    theta = -angulo
    c, s = np.cos(theta), np.sin(theta)
    R = np.array(((c, -s), (s, c)))
    points_rot=[]
    for i in points:
        points_rot.append(np.dot(R,i))    
    xrott=[]
    yrott=[]    
    for i in points_rot: 
        xrott.append(i[0])
        yrott.append(i[1])
 # Traslado en "x" y en "y" 
    yfi=np.array(yrott)-np.mean([yrott[-1],yrott[0]])
    xfi=np.array(xrott)-(max(xrott)+min(xrott))/2
    return xfi, yfi
#%%
xfi, yfi =rotar_trasladar (x_posta4,y_posta4)   
  #%% # Puntos finales sin ningun borde que no quiera 
fig, ax = plt.subplots(figsize=(2*6.4, 2*4.8))

ax.plot(xfi,yfi,"ro",markersize=2)
plt.xlabel(r'$x(\mu m )$', fontsize=18)

plt.ylabel(r'$y(\mu m )$', fontsize=18)    
#%% CHEQUEO SI HAY PUNTOS BIVALUADOS

len(xfi)==len(set(xfi))
# Si da TRUE  es que no hay puntos bivaluados, si da FALSE es que no hay
# Para este perfil dio TRUE. No hay puntos bivaluados.

# Entonces lo que hay que hacer es directo guardar los puntos como coordenadas (x, y)
#%% #guardo los datos resultantes
xarray = np.array(xfi)
yarray = np.array(yfi)

data = np.array([xarray, yarray])
data = data.T

with open("euglenido_paper_listo.txt", 'w+') as datafile_id:
#here you open the ascii file
    np.savetxt(datafile_id, data, fmt=['%.20f','%.20f'])  # 20 es el numero de decimales que graba
#%%
# A continuacion pongo algunas estrategias extras para hacer si HAY puntos bivaluados 

# 1) Quedarse con el promedio de los valores bivaluados 
xf=[]
yf=[]
for i in sorted(set(xfi)):
    indices=np.argwhere(yfi==i)
    ys=yfi[indices]
    xf.append(i)
    yf.append(float(np.mean(ys))) #cambiando np.mean por max o min me quedo con el maximo o minimo
#%%
# 2) Tambien encontre "filtros" muy comunes en analisis de se;ales. No los use, pero sirven para hacer mas suave/ menos picudo el perfil.    
# Probemoslo con el perfil que acabamos de armar
    

def smooth(x,window_len=20,window='bartlett'):
    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    
    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[int(window_len/2-1):-int(window_len/2)]
  

yhat=smooth(yfi)
  #%% # Comparacion entre puntos originales y con filtro smooth
fig, ax = plt.subplots(figsize=(2*6.4, 2*4.8))

ax.plot(xfi,yfi,"ro",markersize=2, label="Sin filtro") #sin filtro
ax.plot(xfi,yhat,"bo",markersize=2, label= "Con filtro Smooth") #con filtro
plt.xlabel(r'$x(\mu m )$', fontsize=18)
plt.legend(loc="best")
plt.ylabel(r'$y(\mu m )$', fontsize=18)  

#%%   3) Otro filtro 
from scipy.signal import savgol_filter
yhat2 = savgol_filter(yfi, 51, 3) # window size 51, polynomial order 3
  #%% # Comparacion entre puntos originales y con filtro savgol
fig, ax = plt.subplots(figsize=(2*6.4, 2*4.8))

ax.plot(xfi,yfi,"ro",markersize=2, label="Sin filtro") #sin filtro
ax.plot(xfi,yhat2,"bo",markersize=2, label= "Con filtro Savgol") #con filtro
plt.xlabel(r'$x(\mu m )$', fontsize=18)
plt.legend(loc="best")
plt.ylabel(r'$y(\mu m )$', fontsize=18)  


#%% Esta funcion es para interpolar los puntos
# Probemos interpolando por ejemplo los puntos  sin filtro
g = interpolate.interp1d(xfi, yfi)
xnew = np.linspace(xfi[0], xfi[-1], 10**4)
ynew = g(xnew)   # use interpolation function returned by `interp1d`
  #%% # Comparacion entre puntos originales y su interpolacion
fig, ax = plt.subplots(figsize=(2*6.4, 2*4.8))

ax.plot(xfi,yfi,"ro",markersize=3, label="Puntos originales") #sin filtro
ax.plot(xnew,ynew,"b",markersize=2, label= "Interpolacion") #con filtro
plt.xlabel(r'$x(\mu m )$', fontsize=18)
plt.legend(loc="best")
plt.ylabel(r'$y(\mu m )$', fontsize=18)  

# Vemos que la interpolacion tiene una parte muy picuda alrededor de x=-0.4 y x=-0.2
#%% Interpolemos con algun filtro a ver si se soluciona 

g = interpolate.interp1d(xfi, yhat)
xnew = np.linspace(xfi[0], xfi[-1], 10**4)
ynew = g(xnew)   # use interpolation function returned by `interp1d`
  #%% # Comparacion entre puntos originales e interpolacion con filtro smooth
fig, ax = plt.subplots(figsize=(2*6.4, 2*4.8))

ax.plot(xfi,yfi,"ro",markersize=3, label="Puntos originales") #sin filtro
ax.plot(xnew,ynew,"b",markersize=2, label= "Interpolacion de puntos con filtro smooth") #con filtro
plt.xlabel(r'$x(\mu m )$', fontsize=18)
plt.legend(loc="best")
plt.ylabel(r'$y(\mu m )$', fontsize=18)  

# Aca vemos que el filtro soluciono las partes "feas" de la interpolacion
#%% # Supongamos ahora que queremos guardar los puntos correspondientes a la celda anterior (la interpolacion de los puntos con filtro smooth, es decir, los puntos en azul de grafico anterior)
# Roto y traslado
xfi2, yfi2 =rotar_trasladar (xnew,ynew) 

#%% #guardo los datos resultantes
xarray = np.array(xfi2)
yarray = np.array(yfi2)

data = np.array([xarray, yarray])
data = data.T

with open("euglenido_paper_listo_con_filtro.txt", 'w+') as datafile_id:
#here you open the ascii file
    np.savetxt(datafile_id, data, fmt=['%.20f','%.20f'])  # 20 es el numero de decimales que graba