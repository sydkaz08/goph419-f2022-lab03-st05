#!/usr/bin/env python
# coding: utf-8

# In[50]:


import numpy as np
import scipy as sp
import sympy as smp
import matplotlib.pyplot as plt
import scipy.integrate as spi
import math


# In[51]:


#Inputs V, returns V
def dzdt(V): #f1
    return V

#inputs V,Z, returns function dvdt
def dvdt(V,Z): #f2
    H=20
    mu=1.827*10**(-5)
    g0=9.811636
    g_p=3.086*10**(-6)
    r=1.5/2*0.01
    rho=7800
    m=(4/3)*np.pi*r**3*rho
    cd=6*np.pi*mu*r/m
    return g0-cd*V+g_p*Z
H=20
mu=1.827*10**(-5)
g0=9.811636
g_p=3.086*10**(-6)
r=1.5/2*0.01
rho=7800
m=(4/3)*np.pi*r**3*rho
cd=6*np.pi*mu*r/m
print(cd)


# In[52]:


def euler(f1, f2, H, dt, z0, v0, t0):
    t=np.zeros(1)
    Z=np.zeros(1)
    V=np.zeros(1)
    Z[0]=z0
    V[0]=v0
    t[0]=t0
    i=0
    while not Z[i]>=H:
        temp_V=V[i]+f2(V[i],Z[i])*dt
        V=np.insert(V,i+1,temp_V)
        
        temp_Z=Z[i]+f1(V[i+1])*dt
        Z=np.insert(Z,i+1,temp_Z)
        
        
        temp_t=(i+1)*dt
        t=np.insert(t,i+1,temp_t)
        
        i=i+1
    e_a=abs((H-max(Z))/H)
    return V, Z, t, e_a


# In[53]:


#Input dzdt, dvdt, height H, step size dt, inital conditions z0, v0, t0
#Returns velocity, height, time, and error
def RungeKutta4(f1, f2, H, dt, z0, v0, t0):
    t=np.zeros(1)
    Z=np.zeros(1)
    V=np.zeros(1)
    Z[0]=z0
    V[0]=v0
    t[0]=t0
    i=0
    while not Z[i]>=H:
        k1=f2(V[i],Z[i])
        l1=f1(V[i])
        k2=f2(V[i]+0.5*k1*dt,Z[i]+0.5*l1*dt)
        l2=f1(V[i]+0.5*k1*dt)
        k3=f2(V[i]+0.5*k2*dt,Z[i]+0.5*l2*dt)
        l3=f1(V[i]+0.5*k2*dt)
        k4=f2(V[i]+k3*dt,Z[i]+l3*dt)
        l4=f1(V[i]+k3*dt)
        
        temp_V=V[i]+(1/6)*(k1+2*k2+2*k3+k4)*dt
        V=np.insert(V,i+1,temp_V)
        ####################################
       
        
        temp_Z=Z[i]+(1/6)*(l1+2*l2+2*l3+l4)*dt
        Z=np.insert(Z,i+1,temp_Z)
        #####################################
        
        temp_t=(i+1)*dt
        t=np.insert(t,i+1,temp_t)
        
        i=i+1
    e_a=abs((H-max(Z))/H)
    return V, Z, t, e_a


# In[54]:


#NOTE COMPARE ALL METHODS TO z=(1/2)(g*t^2), z=height
H=20
dt=0.5 #Note this is h in above equations in markdown
z0=0
t0=0
v0=0

V, Z, t, e_a=euler(dzdt, dvdt, 20, 0.001, z0, t0, v0)
plt.plot(t,Z,label="True Solution",color="black",linewidth="8")
#print(max(Z),max(t))


V, Z, t, e_a=euler(dzdt, dvdt, 20, dt, z0, t0, v0)
plt.plot(t,Z,".-",label="Euler",color="blue")
#print(max(Z),max(t))

V, Z, t, e_a=RungeKutta4(dzdt, dvdt, H, dt, z0, t0, v0)
plt.plot(t,Z,".-",label="4th order Runge-Kutta",color="red")
#print(max(Z),max(t))

plt.legend()
plt.gca().invert_yaxis()
plt.xlabel("Time (s)")
plt.ylabel("Height (m)")
plt.title("Height vs Time for H=20m, step size of 0.2s")
plt.savefig('/users/alikazmi/Desktop/GOPH/goph419-f2022-lab03-st05/source/Figures/Figure(1)')
plt.show()

V, Z, t, e_a=euler(dzdt, dvdt, 10, 0.001, z0, t0, v0)
plt.plot(t,Z,label="True Solution",color="black",linewidth="8")
#print(max(Z),max(t))

V, Z, t, e_a=euler(dzdt, dvdt, 10, dt, z0, t0, v0)
plt.plot(t,Z,".-",label="Euler",color="blue")
#print(max(Z),max(t))

V, Z, t, e_a=RungeKutta4(dzdt, dvdt, 10, dt, z0, t0, v0)
plt.plot(t,Z,".-",label="4th order Runge-Kutta",color="red")
#print(max(Z),max(t))

plt.legend()
plt.gca().invert_yaxis()
plt.xlabel("Time (s)")
plt.ylabel("Height (m)")
plt.title("Height vs Time for H=10m, step size of 0.2s")
plt.savefig('/users/alikazmi/Desktop/GOPH/goph419-f2022-lab03-st05/source/Figures/Figure(2)')
plt.show()

V, Z, t, e_a=euler(dzdt, dvdt, 40, 0.001, z0, t0, v0)
plt.plot(t,Z,label="True Solution",color="black",linewidth="8")
print(max(Z),max(t))

V, Z, t, e_a=euler(dzdt, dvdt, 40, dt, z0, t0, v0)
plt.plot(t,Z,".-",label="Euler",color="blue")
print(max(Z),max(t))


V, Z, t, e_a=RungeKutta4(dzdt, dvdt, 40, dt, z0, t0, v0)
plt.plot(t,Z,".-",label="4th order Runge-Kutta",color="red")
print(max(Z),max(t))

plt.legend()
plt.gca().invert_yaxis()
plt.xlabel("Time (s)")
plt.ylabel("Height (m)")
plt.title("Height vs Time for H=40m, step size of 0.2s")
plt.savefig('/users/alikazmi/Desktop/GOPH/goph419-f2022-lab03-st05/source/Figures/Figure(3)')
plt.show()


# In[55]:


H=20
mu=1.827*10**(-5)
g0=9.811636
g_p=3.086*10**(-6)
r=1.5/2*0.01
rho=7800
m=(4/3)*np.pi*r**3*rho
cd=6*np.pi*mu*r/m

# print(math.sqrt(H/(1/2*g0)))

H=20
dt=0.05 #Note this is h in above equations in markdown
z0=0
t0=0
v0=0
i=0.001
E1=np.array([])
E3=np.array([])
k=0
while i<10.001:
    V, Z, t, e_a1=euler(dzdt, dvdt, H, i, z0, t0, v0)
    V, Z, t, e_a3=RungeKutta4(dzdt, dvdt, H, i, z0, t0, v0)
    E1=np.insert(E1,k,e_a1)
    E3=np.insert(E3,k,e_a3)
    i=i+0.001
    k=k+1
i=np.arange(0.001,i,0.001) #do -0.001
print(len(i),len(E1))
plt.plot(i,E1,label="Euler",color="blue")
plt.plot(i,E3,label="Runge-kutta 4",color="red")
plt.xlabel("Step Size (s)")
plt.ylabel("Error in H (%)")
plt.legend()
plt.title("Change in Error of $H$ with step size")
plt.savefig('/users/alikazmi/Desktop/GOPH/goph419-f2022-lab03-st05/source/Figures/Figure(4)')
plt.show()


# In[56]:


#ALL THIS SI DOING IS ZOOMING IN ON THE AREA FROM 0 to 2 seconds


H=20
mu=1.827*10**(-5)
g0=9.811636
g_p=3.086*10**(-6)
r=1.5/2*0.01
rho=7800
m=(4/3)*np.pi*r**3*rho
cd=6*np.pi*mu*r/m

# print(math.sqrt(H/(1/2*g0)))

H=20
dt=0.05 #Note this is h in above equations in markdown
z0=0
t0=0
v0=0
i=0.001
E1=np.array([])
E3=np.array([])
k=0
while i<4.001:
    V, Z, t, e_a1=euler(dzdt, dvdt, H, i, z0, t0, v0)
    V, Z, t, e_a3=RungeKutta4(dzdt, dvdt, H, i, z0, t0, v0)
    E1=np.insert(E1,k,e_a1)
    E3=np.insert(E3,k,e_a3)
    i=i+0.001
    k=k+1
i=np.arange(0.001,i,0.001) #do -0.001
print(len(i),len(E1))
plt.plot(i,E1,label="Euler",color="blue")
plt.plot(i,E3,label="Runge-kutta 4",color="red")
plt.xlabel("Step Size (s)")
plt.ylabel("Error in H (%)")
plt.legend()
plt.title("Change in Error of $H$ with step size between 0 and 2")
plt.savefig('/users/alikazmi/Desktop/GOPH/goph419-f2022-lab03-st05/source/Figures/Figure(5)')
plt.show()


# In[57]:


#SMALL CHANGE IN g0, g', and cd*

H=20
dt=0.01 #Note this is h in above equations in markdown
z0=0
t0=0
v0=0

def dzdt(V): #f1
    return V

i=0
for i in range(1,11):
    H=20
    mu=1.827*10**(-5)
    g01=9.811636
    g0=g01+g01*0.01*i
    g_p=3.086*10**(-6)
    r=1.5/2*0.01
    rho=7800
    m=(4/3)*np.pi*r**3*rho
    cd=6*np.pi*mu*r/m
    def dvdt(V,Z): #f2
        return g0-cd*V+g_p*Z
    V, Z, t, e_a=euler(dzdt, dvdt, H, dt, z0, t0, v0)
    plt.plot(t,Z,label=f"i={i*1}")
plt.title("Change in $g_0$")
plt.xlabel("Time (s)")
plt.ylabel("Height (m)")
plt.legend()
plt.gca().invert_yaxis()
plt.savefig('/users/alikazmi/Desktop/GOPH/goph419-f2022-lab03-st05/source/Figures/Figure(6)')
plt.show()


i=0
for i in range(1,11):
    H=20
    mu=1.827*10**(-5)
    g0=9.811636
    g_p1=3.086*10**(-6)
    g_p=g_p1+0.01*i*g_p1
    r=1.5/2*0.01
    rho=7800
    m=(4/3)*np.pi*r**3*rho
    cd=6*np.pi*mu*r/m
    def dvdt(V,Z): #f2
        return g0-cd*V+g_p*Z
    V, Z, t, e_a=euler(dzdt, dvdt, H, dt, z0, t0, v0)
    plt.plot(t,Z,label=f"i={i*1}")
plt.title("Change in $g'$")
plt.xlabel("Time (s)")
plt.ylabel("Height (m)")
plt.legend()
plt.gca().invert_yaxis()
plt.savefig('/users/alikazmi/Desktop/GOPH/goph419-f2022-lab03-st05/source/Figures/Figure(7)')
plt.show()

i=0
for i in range(1,11):
    H=20
    mu=1.827*10**(-5)
    g0=9.811636
    g_p=3.086*10**(-6)
    r=1.5/2*0.01
    rho=7800
    m=(4/3)*np.pi*r**3*rho
    cd1=6*np.pi*mu*r/m
    cd=cd1+i*cd1*0.01
    def dvdt(V,Z): #f2
        return g0-cd*V+g_p*Z
    V, Z, t, e_a=euler(dzdt, dvdt, H, dt, z0, t0, v0)
    plt.plot(t,Z,label=f"i={i*1}")
plt.title("Change in $c^*_D$")
plt.xlabel("Time (s)")
plt.ylabel("Height (m)")
plt.legend()
plt.gca().invert_yaxis()
plt.savefig('/users/alikazmi/Desktop/GOPH/goph419-f2022-lab03-st05/source/Figures/Figure(8)')
plt.show()


# In[58]:


#THIS CODE IS SEEING HOW MANY ITERATIONS IT TAKES TO GET A 10ms DIFFERENCE in g0

H=20
dt=0.01 #Note this is h in above equations in markdown
z0=0
t0=0
v0=0
def dvdt(V,Z): #f2
    H=20
    mu=1.827*10**(-5)
    g0=9.811636
    g_p=3.086*10**(-6)
    r=1.5/2*0.01
    rho=7800
    m=(4/3)*np.pi*r**3*rho
    cd=6*np.pi*mu*r/m
    return g0-cd*V+g_p*Z
V, Z, t, e_a=RungeKutta4(dzdt, dvdt, H, dt, z0, t0, v0)
MAX_T=max(t)
plt.plot(t,Z,label="TRUE")

diff=0
i=0
while not diff>=0.01:
    i=i+1
    H=20
    mu=1.827*10**(-5)
    g01=9.811636
    g0=g01+g01*1*0.01*i
    g_p=3.086*10**(-6)
    r=1.5/2*0.01
    rho=7800
    m=(4/3)*np.pi*r**3*rho
    cd=6*np.pi*mu*r/m
    def dvdt(V,Z): #f2
        return g0-cd*V+g_p*Z
    V, Z, t, e_a=euler(dzdt, dvdt, H, dt, z0, t0, v0)
    max_t=max(t)
    diff=abs(MAX_T-max_t)
plt.plot(t,Z,"--",color="red",label=f"iterations={i}")
print("After",i,"iterations")
print(max_t,MAX_T)
print(g0,9.811636)
plt.title("Change in $g_0$ needed for maximum time change of 10ms")
plt.xlabel("Time (s)")
plt.ylabel("Height (m)")
plt.legend()
plt.gca().invert_yaxis()
plt.savefig('/users/alikazmi/Desktop/GOPH/goph419-f2022-lab03-st05/source/Figures/Figure(9)')
plt.show()


# In[64]:


H=20
dt=0.01 #Note this is h in above equations in markdown
z0=0
t0=0
v0=0
def dvdt(V,Z): #f2
    H=20
    mu=1.827*10**(-5)
    g0=9.811636
    g_p=3.086*10**(-6)
    r=1.5/2*0.01
    rho=7800
    m=(4/3)*np.pi*r**3*rho
    cd=6*np.pi*mu*r/m
    return g0-cd*V+g_p*Z
V, Z, t, e_a=RungeKutta4(dzdt, dvdt, H, dt, z0, t0, v0)
MAX_T=max(t)
plt.plot(t,Z,label="TRUE")

diff=0
i=0
while not diff>=0.01:
    i=i+1
    H=20
    mu=1.827*10**(-5)
    g0=9.811636
    g_p1=3.086*10**(-6)
    g_p=g_p1+10*i*g_p1
    r=1.5/2*0.01
    rho=7800
    m=(4/3)*np.pi*r**3*rho
    cd=6*np.pi*mu*r/m
    def dvdt(V,Z): #f2
        return g0-cd*V+g_p*Z
    V, Z, t, e_a=euler(dzdt, dvdt, H, dt, z0, t0, v0)
    max_t=max(t)
    diff=abs(MAX_T-max_t)
print(max_t,MAX_T)
print(g_p,3.086*10**(-6))
plt.plot(t,Z,"--",color="red",label=f"iterations={i*1000}")
print("After",i,"iterations")
plt.title("Change in $g'$ needed for maximum time change of 10ms")
plt.xlabel("Time (s)")
plt.ylabel("Height (m)")
plt.legend()
plt.gca().invert_yaxis()
plt.savefig('/users/alikazmi/Desktop/GOPH/goph419-f2022-lab03-st05/source/Figures/Figure(10)')
plt.show()


# In[63]:


#THIS CODE IS SEEING HOW MANY ITERATIONS IT TAKES TO GET A 10ms DIFFERENCE in cd*

H=20
dt=0.01 #Note this is h in above equations in markdown
z0=0
t0=0
v0=0
def dvdt(V,Z): #f2
    H=20
    mu=1.827*10**(-5)
    g0=9.811636
    g_p=3.086*10**(-6)
    r=1.5/2*0.01
    rho=7800
    m=(4/3)*np.pi*r**3*rho
    cd=6*np.pi*mu*r/m
    return g0-cd*V+g_p*Z
V, Z, t, e_a=euler(dzdt, dvdt, H, dt, z0, t0, v0)
MAX_T=max(t)
plt.plot(t,Z,label="TRUE")

diff=0
i=0
while not diff>=0.01:
    i=i+1
    H=20
    mu=1.827*10**(-5)
    g0=9.811636
    g_p=3.086*10**(-6)
    r=1.5/2*0.01
    rho=7800
    m=(4/3)*np.pi*r**3*rho
    cd1=6*np.pi*mu*r/m
    cd=cd1+i*cd1*0.01
    def dvdt(V,Z): #f2
        return g0-cd*V+g_p*Z
    V, Z, t, e_a=RungeKutta4(dzdt, dvdt, H, dt, z0, t0, v0)
    max_t=max(t)
    diff=abs(MAX_T-max_t)
print(max_t,MAX_T)
print(cd,6*np.pi*mu*r/m)
print(abs(cd-6*np.pi*mu*r/m)/(6*np.pi*mu*r/m))
plt.plot(t,Z,"--",color="red",label=f"iterations={i}")
print("After",i,"iterations")
plt.title("Change in $c_D^*$ needed for maximum time change of 10ms")
plt.xlabel("Time (s)")
plt.ylabel("Height (m)")
plt.legend()
plt.gca().invert_yaxis()
plt.savefig('/users/alikazmi/Desktop/GOPH/goph419-f2022-lab03-st05/source/Figures/Figure(11)')
plt.show()


# In[62]:


#SEES HOW HEIGHT IMPACTS THE ERROR

mu=1.827*10**(-5)
g0=9.811636
g_p=3.086*10**(-6)
r=1.5/2*0.01
rho=7800
m=(4/3)*np.pi*r**3*rho
cd=6*np.pi*mu*r/m

E=np.zeros(9991)
I=np.zeros(9991)
for i in range(10,10001,1):
    V, Z, t, e_a=euler(dzdt, dvdt, i, .2, z0, t0, v0)
    E[i-10]=e_a
    I[i-10]=i
plt.plot(I,E)
plt.xlabel("Height (m)")
plt.ylabel("Error in height (%)")
plt.title("Impact of error on height")
plt.savefig('/users/alikazmi/Desktop/GOPH/goph419-f2022-lab03-st05/source/Figures/Figure(12)')
plt.show()


# In[ ]:




