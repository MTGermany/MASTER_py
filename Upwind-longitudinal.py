#!/usr/bin/env python
# coding: utf-8

SMALL_VAL=1e-7

# In[1]:


L=5000;
dx=50;
Nx=int(L/dx);
dt=0.5;
T=90.5;
Nt=T/dt;
tau=0.001;  # 0.01 for LWR-like, 1-3 for instabilities
timegap=1.2;
rhomax=0.160;
s0=0; # minimum gap; rho_jam = 1/(1/rhomax+s0) 
gamma=1.5;
delta=0; # addtl gradient term

Vf=110/3.6;
leff=1/rhomax;

# simulation control variables

initIndex=2  # {free=0, congested=1, both=2, testFlowInstability=3}
use_upwind=True  #true: upwind; false: leapfrog
# In[2]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Constants

road_vector = np.linspace(0,L,Nx)

# Transition points

rho_crt = 1/ (Vf*timegap + leff) #critical density "at capacity"
rho_jam = 1/(1/rhomax+s0)    # where Ve=0 but not yet bumper2bumper (rhomax)
flow_max = Vf*rho_crt


# Initialize rho initial array

rho_init = np.zeros((Nx,1))


# Initial condition


if(initIndex==0):
    rho_init[0:int(0.25*Nx)] = 0.001  # last element exclusive
    rho_init[int(0.25*Nx):Nx] = 0.015
elif(initIndex==1):
    rho_init[0:int(0.25*Nx)] = 0.040
    rho_init[int(0.25*Nx):Nx] = 0.060
elif(initIndex==2):
    rho_init[0:int(0.15*Nx)] = 0.00
    rho_init[int(0.15*Nx):int(0.50*Nx)] = 0.160
    rho_init[int(0.50*Nx):Nx] = 0.00
else:
    rho_init[0:int(0.50*Nx)] = 0.040
    rho_init[int(0.50*Nx):int(0.55*Nx)] = 0.042
    rho_init[int(0.55*Nx):Nx] = 0.040

    
    
# continuous transition

if(initIndex<2):
    i1=int(0.25*Nx)
    i2=int(0.35*Nx)
    for i in range(i1,i2):
        rho_init[i]=rho_init[i1-1]+(i-i1)/(i2-i1)*(rho_init[i2]-rho_init[i1-1])

# bring ICs within limits

for i in range(0,Nx):
    if rho_init[i] < SMALL_VAL:
        rho_init[i] = SMALL_VAL
    elif rho_init[i] > rhomax:
        rho_init[i] = rhomax

#print('rho_init=',rho_init)




# In[3]:


#import numpy as np
#import matplotlib.pyplot as plt

# Congestion wave speed
w = -flow_max / (rho_jam - rho_crt)

print('w=',w)
print( "rho_crt=",rho_crt)
print( "rho_jam=",rho_jam)

# Initialize an array to hold the flow values

flow_init = np.zeros((Nx,1))

# Calculate the flow for each density
for i in range(0,len(rho_init)):
    if rho_init[i] <= rho_crt:
        flow_init[i] = rho_init[i] * Vf
    elif rho_init[i] < rho_jam:
        flow_init[i] = flow_max + w * (rho_init[i] - rho_crt)
    else:
        flow_init[i] =0

# flow_init[0] = 0

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(1000*rho_init, 3600*flow_init, label='Flow (q) vs Density (k)')
plt.xlabel('Density (k) [vehicles/km]')
plt.ylabel('Flow (q) [vehicles/hour]')
plt.title('Triangular Fundamental Diagram')
plt.legend()
plt.grid(True)
#plt.show()
plt.savefig('FD_range.png', dpi=600)


# In[4]:


# Initialize the equilibrium velocity
V_init = np.zeros((Nx,1))  #Nx times 1 matrix
#V_init = np.zeros(Nx)  #Nx array

V_init=flow_init/rho_init
        
    
# plotting the initial velocity
plt.figure(figsize=(10, 6))
plt.plot(road_vector,V_init *3.6)
plt.xlabel('Road Length (m)')
plt.ylabel('Speed (Km/hr)')


# In[5]:


# Calculate anticipation distance
ant_dist = gamma * V_init * timegap
#print('ant_dist=',ant_dist)
print('ant_dist[0]=',ant_dist[0])

# Calculate k
k = np.floor(ant_dist / dx).astype(int) #Nx times 1 matrix, zeroes generally
#print ('k=',k)
ant_rho = np.zeros((Nx,1))

# Anticipated Density

for i in range(0,Nx):
    #ant_rho[i] = rho_init[i+ k[i]] + ((rho_init[i + k[i] + 1] - rho_init[i + k[i]]) * ((ant_dist[i] / dx) - k[i]))

    # enable anticipation distances>1 cell
    
    i1=min(i+ k[i],Nx-1)
    i2=min(i+ k[i]+1,Nx-1)
    ant_rho[i] = rho_init[i1] + ((rho_init[i2] - rho_init[i1]) * ((ant_dist[i] / dx) - k[i]))

# ant_rho[-1] = rho_init[-1] #superfluous
# print ('ant_rho=',ant_rho)

# plotting the anticipated density
plt.figure(figsize=(10, 6))
plt.plot(road_vector,ant_rho*1000, label = 'Anticipated density')
plt.plot(road_vector,rho_init*1000,label = 'Initial density')
plt.legend(loc= 'upper left')
plt.xlabel('Road Length (m)')
plt.ylabel('Density (Veh/km)')
#plt.show()

# In[6]:


# Initialize the Speed with adaptation density and initial flow
V_star = np.zeros((Nx, 1))

#for i in range(0,Nx):
#    V_star[i] = flow_init[i] / ant_rho[i]
V_star=flow_init/ant_rho

# plotting the V star

plt.plot(road_vector,V_star *3.6, label = 'V_star')
plt.plot(road_vector,V_init *3.6, label = 'Equilibrium Velocity')
plt.legend(loc= 'upper right')
plt.xlabel('Road Length (m)')
plt.ylabel('Speed (Km/hr)')


# In[7]:


# Initialize the derivative of speed with respect to density
dV_drho = np.zeros((Nx, 1))

for i in range(0,Nx):
    if rho_init[i]<=rho_crt:
        dV_drho[i] = 0
    else:
        #dV_drho[i] = - flow_max/(rho_init[i]**2)
        dV_drho[i] = - 1/(rho_init[i]**2*timegap)


# In[8]:


# Define display times and initialize related variables

display_time = [30, 60,90] # must be integer for selection

z = 0
count = len(display_time)

# dataframe to store the densities at different time intervals
df = pd.DataFrame()

# Defining the final density , flow  and anticipation density
rho_new = np.zeros((Nx,1))
Q_new = np.zeros((Nx,1))
ant_rho_new = np.zeros((Nx,1))



z = 0
for i in range(0,Nx):
    df.at[i,'road_vector'] = z
    z = z + dx


# In[9]:


df['Initial'] = rho_init


# In[10]:


# Initializing the variables for the simulation

rho_old = rho_init.copy()
u_old = V_init.copy()
Q_old = flow_init.copy()
Ve_star = V_star.copy()
dv = dV_drho.copy()

#########################################################
# The actual simulation
#########################################################

import math
for m in range(0,int(Nt)):
    U1 = Q_old/rho_old
    U2 = (Q_old **2) / rho_old

    for i in range(1, Nx-1):

        if(use_upwind):
            rho_new[i] = rho_old[i] - dt/dx * (Q_old[i] - Q_old[i-1])
            Q_new[i] = Q_old[i] - dt/dx * (U2[i] - U2[i-1]) + (rho_old[i] * Ve_star[i] - Q_old[i])*(1-math.exp(-dt/tau)) - delta*(rho_old[i]**2)*dv[i]* (dt/dx)*(U1[i] - U1[i-1])   # reformulated relax. term allows for lower tau

        else:
            rho_new[i] = rho_old[i] - 0.5*dt/ dx * (Q_old[i+1] - Q_old[i-1])
            Q_new[i] = Q_old[i] - 0.5*dt/dx * (U2[i+1] - U2[i-1]) + (rho_old[i] * Ve_star[i] - Q_old[i])*(1-math.exp(-dt/tau)) - delta*(rho_old[i]**2)*dv[i]* (dt/dx)*(U1[i] - U1[i-1])   # reformulated relax. term allows for lower tau


        # Obviously, Godunov does not work so simple...
        # di=0
        # if (i<Nx-1) and ((rho_old[i-1]+rho_old[i]+rho_old[i+1])/3.>rho_crt):
        #    di=1
        
        # rho_new[i] = rho_old[i] - ((dt / dx) * (Q_old[i+di] - Q_old[i+di-1]))
        # Q_new[i] = Q_old[i] - ((dt / dx) * (U2[i+di] - U2[i+di-1])) + dt*(rho_old[i] * Ve_star[i] - Q_old[i])/tau - delta*(rho_old[i]**2)*dv[i]* (dt/dx)*(U1[i] - U1[i-1])        

    # Von-Neumann Boundary Conditions
    
    rho_new[0] = rho_new[1]                       # rho_new(0) = rho_new(1)
    rho_new[-1] = rho_new[-2]                     # rho_new(end) = rho_new(end-1)
    Q_new[0] = Q_new[1]                           # Q_new(0) = Q_new(1)
    Q_new[-1] = Q_new[-2]                         # Q_new(end) = Q_new(end-1)
    
    
    # #Updated Speed with new density and Flow

    V_update = Q_new / rho_new

            
    # Updated anticipated distance

    #ant_dist_update = gamma * V_update * timegap
    ant_dist_update=np.zeros((Nx,1))
    for i in range(0,Nx):
        ant_dist_update[i]=gamma * Vf * timegap


    # finding the k = floor(anticipated distance / dx )

    k = np.floor(ant_dist_update / dx).astype(int)

    
    # Update the anticipated density
    # enable anticipation distances>1 cell

    for i in range(0,Nx):
        i1=min(i+ k[i],Nx-1)
        i2=min(i+ k[i]+1,Nx-1)
        ant_rho_new[i] = rho_new[i1] + ((rho_new[i2] - rho_new[i1]) * ((ant_dist_update[i] / dx) - k[i]))

    
    #ant_rho_new[-1] = rho_new[-1]       #superfluous
    
    # updating the V_star with new flow and new density
    Ve_star_update = np.zeros((Nx,1))
    for i in range(0,Nx):
        if ant_rho_new[i]<=rho_crt:
            Ve_star_update[i] = Vf
        elif(ant_rho_new[i]<=rho_jam):
            #Ve_star_update[i] = Q_new[i] / ant_rho_new[i] # wrong
            Ve_star_update[i] = 1./(ant_rho_new[i]*timegap) - 1./(timegap*rho_jam)
        else:
            Ve_star_update[i]=0

    
    #updating the derivative of speed with respect to density
    dv_update = np.zeros((Nx,1))
    for i in range(0,Nx):
        if rho_new[i]<=rho_crt:
            dv_update[i] = 0
        else:
            #dv_update[i] = -flow_max / (rho_new[i]**2) #wrong
            dv_update[i] = -1 / (rho_new[i]**2*timegap)
            
    # Defining the limits of the density between 0 and rhomax

    for i in range(0,Nx):
        if rho_new[i] < SMALL_VAL:
            rho_new[i] = SMALL_VAL
        elif rho_new[i] > rhomax:
            rho_new[i] = rhomax
        #else:
        #    rho_new[i] = rho_new[i]

    Q_old = Q_new.copy()
    rho_old = rho_new.copy()
    Ve_star = Ve_star_update.copy()
    dv = dv_update.copy()

    #debug
    
    deb_flag=False
    for i in range(0, Nx-1):
        if (rho_old[i]==0):
            deb_flag=True
    if(deb_flag):
        for i in range(0, Nx-1):
            print("i=",i," m=it=",m,"rho_old[i]=",rho_old[i])
   
    #storing density in a dataframe
    
    if int(m*dt) in display_time: # corrected error
        df[f'{int(m*dt)} sec'] = rho_new.tolist()  # Adjust indexing for storing_den
    


# In[11]:


df = df.apply(lambda col: col.map(lambda x: x[0] if isinstance(x, (np.ndarray, list)) else x))


# In[12]:


plt.figure(figsize=(10, 6))
plt.plot(df['road_vector'].values,df['Initial'].values*1000, label = 'initial')
plt.plot(df['road_vector'].values,df['30 sec'].values*1000, label = '30 sec')
plt.plot(df['road_vector'].values,df['60 sec'].values*1000, label = '60 sec')
plt.plot(df['road_vector'].values,df['90 sec'].values*1000, label = '90 sec')
plt.legend(loc= 'upper left')
plt.xlabel('Road Length (m)')
plt.ylabel('Density (Veh/km)')
plt.grid(True)
plt.savefig('Density_Distribution.png', dpi=600)


# In[ ]:





# In[ ]:





# In[ ]:




