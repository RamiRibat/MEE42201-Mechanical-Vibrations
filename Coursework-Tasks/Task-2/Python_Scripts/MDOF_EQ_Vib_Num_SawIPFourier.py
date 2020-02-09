# Mech Vib Tsak 2 (2DOF, SAW Earthquake Vibration, Numeric Solution) #

#************************************************************************************************

    ## The required Libraries to be imported:

import numpy as np
from numpy import rad2deg
from numpy.linalg import inv, eigh
import matplotlib.pyplot as plt
from math import pi, sqrt, sin, cos, atan, exp
from sympy import symbols
from sympy.physics.vector import init_vprinting, vlatex
init_vprinting(use_latex='mathjax', pretty_print=False)

#************************************************************************************************

    ## Functions
def Y(t):
    Y = np.zeros((2*dof,1))
    dyn = 0
    yn = 0#-Yo/2
    for k in range(1,n):
        dyn_new = dyn - Yo * k*w_eq * cos(k*w_eq*t)/(k*pi);
        dyn = dyn_new
    for m in range(1,n):
        yn_new = yn - Yo * sin(m*w_eq*t)/(m*pi);
        yn = yn_new
    Y[0] = k1*yn + c1*dyn
    return Y

def Yn(t):  
    yn = 0#-Yo/2
    for m in range(1,n):
        yn_new = yn - Yo * sin(m*w_eq*t)/(m*pi);
        yn = yn_new
    Yn = yn
    return Yn

def G(x,t):
    return A_inv.dot( Y(t) - B.dot(x) )

def RK4_step(x,t,dt): #Runge-Kutta
    k1 = G(x, t)
    k2 = G(x+0.5*k1*dt, t+0.5*dt)
    k3 = G(x+0.5*k2*dt, t+0.5*dt)
    k4 = G(x+k3*dt, t+dt)
    #return dt * G(y,t)
    return dt * (k1 + 2*k2 + 2*k3 + k4)/6

#************************************************************************************************

    ## System Parameters:

        ### Materials:
        ### Materials:
E_s = 2.1e11 # Modulus of rigidity of steel [Pa]
E_c = 0.17e11 # Modulus of rigidity of concrete [Pa]
E = 0.2*E_s + 0.8*E_c # Total modulus of rigidity [Pa]
rho = 1000 # RC density [kg/m^3]

        ### Dimensions:
L = 5 # Floor length [m]
W = 5 # Floor width [m]
h = 0.5 # Floor thickness [m]
H  = 5 # Story hight [m]
l = 0.65 # Column length [m]
w = 0.65 # Column width [m]
V = L * W * h # Floor volume

        ### Dynamic Parametrs:
m = rho * V # Floor mass [kg]

I = (l * w**4) / 12 # Column inertia [m^4]
ko = (3 * E * I) / H**3 # Column stiffness (cantilever beam)
k = 4*ko
k1 =k

zeta = 0.05 # Daming ratio of the building
wn = sqrt(k/m) # Natural frequency
wd = sqrt(1-zeta**2)*wn # Damped Natural frequency

cc =  2*sqrt(m*k) # Critical damping
c = zeta * cc # Damping coefficient
c1 = c

#************************************************************************************************

dof = 10

M = np.zeros((dof,dof))
C = np.zeros((dof,dof))
K = np.zeros((dof,dof))
I = np.identity(dof)

for n in range(0,dof):
    M[n,n] = m

for n in range(0,dof-1):
    C[n,n] = 2*c
    C[n,n+1] = -c
    C[n+1,n] = C[n,n+1]
C[dof-1,dof-1] = c 

for n in range(0,dof-1):
    K[n,n] = 2*k
    K[n,n+1] = -k
    K[n+1,n] = K[n,n+1]
K[dof-1,dof-1] = k


A = np.zeros((2*dof,2*dof))
B = np.zeros((2*dof,2*dof))

A[0:dof,0:dof] = M
A[dof:2*dof,dof:2*dof] = I

B[0:dof,0:dof] = C
B[0:dof,dof:2*dof] = K
B[dof:2*dof,0:dof] = -I

A_inv = inv(A)

#************************************************************************************************

    ## Input Data (SAW wave):

n = 5; # number of Fourier terms
freq_eq = 10 # Excitation frequency [Hz]
w_eq = 2 * pi * freq_eq # Cicular frequency [rad/s]
Yo = 5 # Earth quick amplitude [m]

tau = 2*pi/w_eq;
#sim_time = 12

#************************************************************************************************

    ## Simulation:

dt = 1/(32*freq_eq) # Time step
end_time = 2 # Finish time
time = np.arange(0, end_time, dt) # time range

        ### History

Xeq = np.zeros((2*dof,1)) # initial conditions
Xn = np.zeros((len(time),dof)) # Disp History
Vn = np.zeros((len(time),dof)) # Velo History
Seismicy = []

V1 = [] # v response history
V2 = [] # v response history
X1 = [] # x response history
X2 = [] # x response history

        ### Time Response:

#--------------------------------------------------
# numerically integrate the EOMs
for k in range(0,len(time)):
    
    t = time[k]
    X_new = Xeq + RK4_step(Xeq,t,dt)
    Xeq = X_new

    X1.append(Xeq[2]) # Displacemant output history
    X2.append(Xeq[3]) # Displacemant output history
    
    for i in range(0,dof):
        Xn[k,i] = Xeq[dof+i]
    
    Seismicy.append(Yn(t))

    
#************************************************************************************************

    ## Results and Responses:

        ### Abs response x_eq(t):
legend_list = []
#j_list = [1, 2, 5, 10, 25, 50, 100]
j_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
for j in range(0,len(j_list)):
    i = j_list[j]
    plt.plot(time,Xn[0:,i-1])
    leg = 'X'+str(i)
    legend_list.append(leg)

plt.plot(time,Seismicy) 
plt.xlabel('time (s)', fontsize=20)
plt.ylabel('displacement (m)', fontsize=20)
plt.title('Seismic Response of 10DOF Building', fontsize=25)
#plt.legend([legend_list[i] for i in range(0,len(legend_list))], loc='lower right', fontsize=13)
plt.legend([j_list[i] for i in range(0,len(j_list))], loc='lower right', fontsize=13)
plt.grid(True)
plt.show()


plt.show()
