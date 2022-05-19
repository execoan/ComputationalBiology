# -*- coding: utf-8 -*-
"""
Created on Wed May 18 10:03:45 2022

@author: yunus
"""

import numpy as np
import matplotlib.pyplot as plt


s = float(input('Please give nutrient concentration: '))         #nutrient concentration

r_max = 2.1           #maximal per_cell growth rate 

gamma = 0.01                 #number of units nutrient that are consumed to produce one bacteria (divided volume)

K_s  = 0.1           #half velocity constant the value of s r/r_max = 0.5 

initial_cell = float(input("Please give initial bacteria number:"))

print('Inıtıal Cell:', initial_cell)


#10dt = 0.1
t_start = 0             # starttime
t_end = 10        # endtime
#n_steps = int(round((t_end-t_start)/dt))  

def get_monod_gr(s,r_max,K_s):
    rate = r_max*s/(K_s+s)
    return rate

def get_N(monod, N):
    return monod*N

def get_S(g,monod,N):
    return -g*monod*N

from scipy.integrate import solve_ivp

gr = []
tt = []


s0 = 1
def fun(t,Y):
    dY = np.zeros_like(Y)
    
    gr.append(get_monod_gr(Y[0], r_max, K_s))
    tt.append(t)
    
    dY[0] = -gamma * get_monod_gr(Y[0], r_max, K_s) * Y[1] #- D* Y[1] 
    dY[1] =  get_monod_gr(Y[0], r_max, K_s) * Y[1] #+ s0*D - Y[0]*D
    return dY

Y0 = np.array([s,initial_cell])

sol = solve_ivp(fun, [t_start, t_end], Y0, dense_output=True)

time = np.linspace(t_start,t_end,1000)
Y = sol.sol(time)

maks_pop = np.max(Y[1])
r_max_logistic = 2.2

solution_logistic_growth_model = (initial_cell * np.exp(r_max_logistic * time)) / (1+ ((initial_cell/maks_pop)*(np.exp(r_max_logistic * time)-1)))  

plt.plot(time,solution_logistic_growth_model,'-',color='green',label='Lojistik Model')
plt.plot(time,Y[0],'-',label='Besin Yoğunluğu',color='orange')
plt.plot(time,Y[1],'-', color='red',label='Monod Fonksiyonu')
plt.xlabel('Süre[h]')
plt.ylabel('N(t)')
plt.legend()
plt.show()

print(np.min(Y[0]))

print(np.shape(gr))
plt.plot(tt, gr,'.')
plt.ylabel('Büyüme Oranı')
plt.xlabel('Süre(h)')
plt.show()


gamma2 = 1 #grafiğe uyması için deneyde alınan değerleri aldım
K_s2 = 0.1
s02 = 1
r_max2 = 2
D = 1.81

def fun2(t,Y):
    
    dY = np.zeros_like(Y)
    
    gr.append(get_monod_gr(Y[0], r_max2, K_s2))
    tt.append(t)
    
    dY[0] = -gamma2 * get_monod_gr(Y[0], r_max2, K_s2) * Y[1] - D* Y[1] 
    dY[1] =  get_monod_gr(Y[0], r_max2, K_s2) * Y[1] + s02*D - Y[0]*D
    return dY

Y0 = np.array([s,initial_cell])

sol = solve_ivp(fun2, [t_start, t_end], Y0, dense_output=True)

time = np.linspace(t_start,t_end,1000)
Y = sol.sol(time)

maks_pop = np.max(Y[1])


plt.plot(time,Y[0],'-',label='Besin Yoğunluğu',color='orange')
plt.plot(time,Y[1],'-', color='red',label='Monod Fonksiyonu')
plt.xlabel('Süre[h]')
plt.ylabel('Kemostat için N(t)')
plt.legend()
plt.show()

expo = initial_cell * np.exp(r_max_logistic * time)
fig, (ax1, ax2) = plt.subplots(2)

fig.suptitle('Üstel ve Lojistik Büyüme Farkı')
ax1.plot(time, expo)
ax2.plot(time,solution_logistic_growth_model)
#ax1.ticklabel_format(useOffset=False, style='plain')
plt.ylabel('                                       Büyüme Oranı')
plt.xlabel('Süre(h)')
plt.show()
