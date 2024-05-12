"""
Este script arma el sistema lineal para resolver el campo reflejado para 
onda plana incidente. El perfil puede ser ingresado con una funcion analitica
o con coordenadas de una imagen
"""

#%%   grilla de alfas
delta=1/(a*3)
xfinn=2*k_0*(epsilon_1*mu_1)**0.5
ene=int((xfinn+alfa_0)/delta)
alfaspos=[]
for i in range(ene+1):
    alfaspos.append(alfa_0+i*delta)
alfaspos=np.array(alfaspos)
alfasneg=[]
for i in reversed(range(1,ene+1)):
    alfasneg.append(alfa_0-i*delta)
alfasneg=np.array(alfasneg) 
alfas_R=np.concatenate(( alfasneg,alfaspos), axis=None)

#%%  solucion particular de Ax=b
indice=np.argwhere(alfas_R==alfa_0)
start_time = time.time()
br=np.zeros((len(alfas_R)),dtype=complex)
br[np.r_[0:indice, indice+1:len(alfas_R)]] = np.array(list(map(b_R,alfas_R[np.r_[0:indice, indice+1:len(alfas_R)]])))
br[indice]=2*np.pi*(sigma_2*beta_alfa(alfa_0,epsilon_1,mu_1)/sigma_1-beta_alfa(alfa_0,epsilon_2,mu_2))*(D(0,beta_alfa(alfa_0,epsilon_2,mu_2)-beta_alfa(alfa_0,epsilon_1,mu_1))-D(0,beta_alfa(alfa_0,epsilon_2,mu_2)+beta_alfa(alfa_0,epsilon_1,mu_1)))
print("--- %s seconds ---" % (time.time() - start_time))

#%% matriz
m_R=np.zeros((len(alfas_R),len(alfas_R)),dtype=complex)
start_time = time.time()
for i in range(len(alfas_R)):
    m_R[i,np.r_[0:i, i+1:len(alfas_R)]]=delta*v_R1(alfas_R[i],alfas_R[np.r_[0:i, i+1:len(alfas_R)]])
    m_R[i,i]=(delta*D(0,beta_alfa(alfas_R[i],epsilon_2,mu_2)-beta_alfa(alfas_R[i],epsilon_1,mu_1))-delta*a/(2*np.pi)+1)*(-beta_alfa(alfas_R[i],epsilon_2,mu_2)-beta_alfa(alfas_R[i],epsilon_1,mu_1)*sigma_2/sigma_1)
    if i%100 == 0:
        t_c = int((time.time() - start_time)/60)
        print('Ya calculamos {} filas y llevamos {} minutos de cálculo'.format(i,t_c))

print("--- %s seconds ---" % (time.time() - start_time))

#%% resuelvo el sistema lineal
x_R=np.linalg.solve(m_R,br)
np.allclose(np.dot(m_R, x_R), br)
#%% guardo los datos
np.savez(r"nombre_del_archivo",matriz=m_R,particular=br,solucion=x_R,alfas_R=alfas_R) 

#%% |R|^2
eje_y_2=[]
for i in x_R:
    eje_y_2.append(np.absolute(i)**2)

#%%
alfa_k1=alfas_R/(k_0*(epsilon_1*mu_1)**0.5)

#%%  curva de |R|^2 vs alfa/k_1
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

fig, ax = plt.subplots(figsize=(12.8,9.6))

ax.plot(alfa_k1,np.array(eje_y_2),'bo',markersize=3)
plt.xlabel(r'$\alpha/k_{1}  $', fontsize=32, labelpad=20)
plt.ylabel(r'$\left | \tilde{R}(\alpha ) \right |^{2}$', fontsize=32, labelpad=20)
#plt.xticks(np.concatenate((np.arange(-2, 2, step=1), 2), axis=None), fontsize=30) esto hay que setearlo acorde al grafico en particular
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(30) 
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(30)    
ax.tick_params(which='both', width=2,length=7, color='k')    
#ax.xaxis.set_minor_locator(MultipleLocator(0.5)) #setear de acuerdo al grafico
#ax.yaxis.set_minor_locator(MultipleLocator(0.0001)) #setear de acuerdo al grafico
ax.annotate(r'pol. s',(0.06,0.85),xycoords='axes fraction',fontsize=32)
plt.tight_layout() 
#%% #para guardar la figura
fig.savefig('R.png', format='png', dpi=500)

#%% DISTRIBUCION ANGULAR  DE POTENCIA
from itertools import repeat
betas_R=np.array(list(map(beta_alfa,alfas_R,repeat(epsilon_1),repeat(mu_1))))

R_abs_cuad=np.array(eje_y_2)
dp=betas_R*R_abs_cuad
dp=dp/(beta_alfa(alfa_0,epsilon_1,mu_1)*2*np.pi)

#%% curva de dp/dalfa vs alfa/k_1
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
fig, ax = plt.subplots(figsize=(12.8,9.6))

ax.plot(alfa_k1,dp,'bo',markersize=3)
ax.yaxis.set_major_locator(plt.MaxNLocator(6))
ax.xaxis.set_major_locator(plt.MaxNLocator(7))
plt.xlabel(r'$\sin\hspace{0.05} \theta_{s1} $', fontsize=32, labelpad=20)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(30) 
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(30) 
ax.tick_params(which='both', width=2,length=7, color='k')    
plt.ylabel(r'$dP^1/d\alpha $', fontsize=32, labelpad=20)

plt.xlim(-1,1) #esto es necesario, asi se selecciona 
#ax.xaxis.set_minor_locator(MultipleLocator(0.15))  # hay que setear de acuerdo al grafico en particular
#ax.yaxis.set_minor_locator(MultipleLocator(0.125))  # hay que setear de acuerdo al grafico en particular 
ax.annotate(r'pol. s',(0.06,0.85),xycoords='axes fraction',fontsize=32)
plt.tight_layout()

#%%



