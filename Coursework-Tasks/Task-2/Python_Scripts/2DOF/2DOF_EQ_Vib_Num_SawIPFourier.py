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
    Y = np.zeros((4,1))
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
m1 = rho * V # Floor mass [kg]
m2 = m1

I = (l * w**4) / 12 # Column inertia [m^4]
ko = (3 * E * I) / H**3 # Column stiffness (cantilever beam)
k1 = 4*ko
k2 = k1

zeta = 0.05 # Daming ratio of the building
wn = sqrt(k1/m1) # Natural frequency
wd = sqrt(1-zeta**2)*wn # Damped Natural frequency

cc =  2*sqrt(m1*k1) # Critical damping
c1 = zeta * cc # Damping coefficient
c2 = c1

M = np.matrix([[m1,0],
               [0,m2]])
C = np.matrix([[c1+c2,-c2],
              [-c2,c2]])
K = np.matrix([[k1+k2,-k2],
              [-k2,k2]])

A= np.zeros((4,4))
A[0:2,0:2] = M
A[2:4,2:4] = np.eye(2)

B= np.zeros((4,4))
B[0:2,0:2] = C
B[0:2,2:4] = K
B[2:4,0:2] = -np.eye(2)

R = np.zeros((4,4))
R[0:2,0:2] = C
B[0:2,2:4] = K

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

        ### Initial Conditions:
v1o = 0.0
v2o = 0.0
x1o = 0.0
x2o = 0.0
x = np.array([[v1o],[v2o],[x1o],[x2o]])

        ### History
V1 = [] # v response history
V2 = [] # v response history
X1 = [] # x response history
X2 = [] # x response history
Yn_history = []
Seismicy = []
Seismicdy = []
        ### Time Response:
#--------------------------------------------------
#Yn_history= [] # SAW input
#for t in time:
#    yn = Yo/2
 #   for m in range(1,n):
  #      yn = yn - Yo * sin(m*w_eq*t)/(m*pi);
   # yn
    #Yn_history.append(yn);
#--------------------------------------------------
for t in time:
    x_new = x + RK4_step(x,t,dt) #Total response, vector & Matrices
    x = x_new

    V1.append(x[0]) # Velocity output history
    V2.append(x[1]) # Velocity output history
    X1.append(x[2]) # Displacemant output history
    X2.append(x[3]) # Displacemant output history
    Seismicy.append(Yn(t))
    #Seismicy.append(Y(t)[2])

    
#************************************************************************************************

    ## Results and Responses:

        ### Abs response x_eq(t):
fig1 = plt.figure()
rect = fig1.patch

graphEQ = fig1.add_subplot(1,1,1)
graphEQ.plot(time, X1, linewidth=2)
graphEQ.plot(time, X2, linewidth=2)
#graphEQ.plot(time, V1, linewidth=2)
#graphEQ.plot(time, V2, linewidth=2)
#graphEQ.plot(time, Seismicy, linewidth=1)
graphEQ.set_title(' "SAW Seismic System Response" ', fontsize=25)
graphEQ.set_xlabel('Time [s]', fontsize=20)
graphEQ.legend(['x1(t) [m]', 'x2(t) [m]', 'y(t) [m]'], fontsize=13, loc='lower right')
graphEQ.grid(True)



plt.show()
