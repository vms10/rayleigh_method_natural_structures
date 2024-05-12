"""
Este script arma el sistema lineal para resolver el campo transmitido para 
haz gaussiano incidente. El perfil puede ser ingresado con una funcion analitica
o con coordenadas de una imagen
"""

#%%   grilla de alfas
delta=1/(a*3)
xfinn=1.5*k_0*(epsilon_2*mu_2)**0.5
ene=int((xfinn+alfa_0)/delta)
alfaspos=[]
for i in range(ene+1):
    alfaspos.append(alfa_0+i*delta)
alfaspos=np.array(alfaspos)
alfasneg=[]
for i in reversed(range(1,ene+1)):
    alfasneg.append(alfa_0-i*delta)
alfasneg=np.array(alfasneg) 
alfas=np.concatenate(( alfasneg,alfaspos), axis=None)

#%%
start_time = time.time()
bt=np.zeros((len(alfas)),dtype=complex)
bt[:] = np.array(list(map(b_T,alfas)))
print("--- %s seconds ---" % (time.time() - start_time))

#%% matriz
m_T=np.zeros((len(alfas),len(alfas)),dtype=complex)
t_0 = time.time()
for i in range(len(alfas)):
    m_T[i,np.r_[0:i, i+1:len(alfas)]]=delta*v_T1(alfas[i],alfas[np.r_[0:i, i+1:len(alfas)]])
    m_T[i,i]=(delta*D(0,beta_alfa(alfas[i],epsilon_2,mu_2)-beta_alfa(alfas[i],epsilon_1,mu_1))-delta*a/(2*np.pi)+1)*(-beta_alfa(alfas[i],epsilon_2,mu_2)-beta_alfa(alfas[i],epsilon_1,mu_1)*sigma_2/sigma_1)    
    if i%(100) == 0:
        t_c = int((time.time() - t_0)/60)
        print('Ya calculamos {} filas y llevamos {} minutos de cálculo'.format(i,t_c))

print("--- %s seconds ---" % (time.time() - t_0))

#%%   resuelvo el sistema lineal
x_T=np.linalg.solve(m_T,bt)
np.allclose(np.dot(m_T, x_T), bt)
#%%  guardo los datos
np.savez(r"nombre_del_archivo",matriz=m_T,particular=bt,solucion=x_T,alfas=alfas) 

#%%
eje_y=[]
for i in x_T:
    eje_y.append(np.absolute(i)**2)

#%%
alfa_k2=alfas/(k_0*(epsilon_2*mu_2)**0.5)

#%% Distribucion RESTANDO lo especular para cada A(alfa)
T = lambda alpha: (2*sigma_2*beta_alfa(alpha,epsilon_1,mu_1)/sigma_1)/(sigma_2*beta_alfa(alpha,epsilon_1,mu_1)/sigma_1+beta_alfa(alpha,epsilon_2,mu_2)) 

resta=[]
for i in range(len(x_T)):
    resta.append(np.absolute(x_T[i]-T(alfas[i])*A_alfa(alfas[i]))**2)
resta=np.array(resta)

from itertools import repeat
betas_T=np.array(list(map(beta_alfa,alfas,repeat(epsilon_2),repeat(mu_2))))
put=resta*betas_T
put=sigma_1*put/(beta_alfa(alfa_0,epsilon_1,mu_1)*sigma_2)
#%% curva de intensidad rugosidades vs alfa/k_2
# Es la que mas informacion contiene
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
fig, ax = plt.subplots(figsize=(12.8,9.6))

ax.plot(alfa_k2,put,'bo',markersize=3)
ax.yaxis.set_major_locator(plt.MaxNLocator(6))
ax.xaxis.set_major_locator(plt.MaxNLocator(7))
plt.xlabel(r'$\sin\hspace{0.05} \theta_{s2} $', fontsize=32, labelpad=20)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(30) 
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(30) 
ax.tick_params(which='both', width=2,length=7, color='k')    
plt.ylabel(r'$Intensidad rugosidades$', fontsize=32, labelpad=20)

plt.xlim(-1,1) #esto es necesario, asi se selecciona 
#ax.xaxis.set_minor_locator(MultipleLocator(0.15))  # hay que setear de acuerdo al grafico en particular
#ax.yaxis.set_minor_locator(MultipleLocator(0.125))  # hay que setear de acuerdo al grafico en particular 
ax.annotate(r'pol. s',(0.06,0.85),xycoords='axes fraction',fontsize=32)
plt.tight_layout()

#%%  curva de |T|^2 vs alfa/k_2
# no tiene mucha info 
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

fig, ax = plt.subplots(figsize=(12.8,9.6))

ax.plot(alfa_k2,np.array(eje_y),'bo',markersize=3)
plt.xlabel(r'$\alpha/k_{2}  $', fontsize=32, labelpad=20)
plt.ylabel(r'$\left | \tilde{T}(\alpha ) \right |^{2}$', fontsize=32, labelpad=20)
#plt.xticks(np.concatenate((np.arange(-2, 2, step=1), 2), axis=None), fontsize=30) esto hay que setearlo acorde al grafico en particular
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(30) 
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(30)    
ax.tick_params(which='both', width=2,length=7, color='k')    
#ax.xaxis.set_minor_locator(MultipleLocator(0.5)) #setear de acuerdo al grafico
#ax.yaxis.set_minor_locator(MultipleLocator(0.0001)) #setear de acuerdo al grafico
ax.annotate(r'pol. s',(0.06,0.85),xycoords='axes fraction',fontsize=32)
plt.ylim()
plt.tight_layout() 
#%% #para guardar la figura
fig.savefig('T.png', format='png', dpi=500)

#%% DISTRIBUCION ANGULAR  DE POTENCIA

from itertools import repeat
betas_T=np.array(list(map(beta_alfa,alfas,repeat(epsilon_2),repeat(mu_2))))

T_abs_cuad=np.array(eje_y)
dp=betas_T*T_abs_cuad
dp=sigma_1*dp/(beta_alfa(alfa_0,epsilon_1,mu_1)*sigma_2*2*np.pi)

#%% curva de dp/dalfa vs alfa/k_2
# no tiene mucha info 
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
fig, ax = plt.subplots(figsize=(12.8,9.6))

ax.plot(alfa_k2,dp,'bo',markersize=3)
ax.yaxis.set_major_locator(plt.MaxNLocator(6))
ax.xaxis.set_major_locator(plt.MaxNLocator(7))
plt.xlabel(r'$\sin\hspace{0.05} \theta_{s2} $', fontsize=32, labelpad=20)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(30) 
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(30) 
ax.tick_params(which='both', width=2,length=7, color='k')    
plt.ylabel(r'$dP^2/d\alpha $', fontsize=32, labelpad=20)

plt.xlim(-1,1) #esto es necesario, asi se selecciona 
#ax.xaxis.set_minor_locator(MultipleLocator(0.15))  # hay que setear de acuerdo al grafico en particular
#ax.yaxis.set_minor_locator(MultipleLocator(0.125))  # hay que setear de acuerdo al grafico en particular 
ax.annotate(r'pol. s',(0.06,0.85),xycoords='axes fraction',fontsize=32)
plt.tight_layout()
