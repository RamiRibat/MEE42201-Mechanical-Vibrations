# Mech Vib Tsak 2 (SDOF, SAW Earthquake Vibration, Numeric Solution) #

#************************************************************************************************

    ## The required Libraries to be imported:

import numpy as np
from numpy import rad2deg
from numpy.linalg import inv
import matplotlib.pyplot as plt
from math import pi, sqrt, sin, cos, atan, exp
from sympy import symbols
from sympy.physics.vector import init_vprinting, vlatex
init_vprinting(use_latex='mathjax', pretty_print=False)

#************************************************************************************************

    ## Functions
def Y(t):
    Y = np.zeros((2,1))
    dyn = 0
    yn = 0#-Yo/2
    for k in range(1,n):
        dyn_new = dyn - Yo * k*w_eq * cos(k*w_eq*t)/(k*pi);
        dyn = dyn_new
    for m in range(1,n):
        yn_new = yn - Yo * sin(m*w_eq*t)/(m*pi);
        yn = yn_new
    Y[0] = k*yn + c*dyn
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
m = rho * V # Floor mass [kg]
I = (l * w**4) / 12 # Column inertia [m^4]
ko = (3 * E * I) / H**3 # Column stiffness (cantilever beam)
k = 4*ko # Equivelant columns stiffnesses
zeta = 0.05 # Daming ratio of the building
wn = sqrt(k/m) # Natural frequency
wd = sqrt(1-zeta**2)*wn # Damped Natural frequency
cc = 2*sqrt(m*k) # Critical damping
c = zeta * cc # Equivelant damping coefficient

A = np.array([[m,0],[0,1]])
B = np.array([[c,k],[-1,0]])
R = np.array([[c,k],[0,0]])
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
vo = 0.0
xo = 0.0
x = np.array([[vo], [xo]])

        ### History
V = [] # v response history
X = [] # x response history
SeismicY = []

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

    V.append(x[0]) # Velocity output history
    X.append(x[1]) # Displacemant output history
    SeismicY.append(Yn(t))

    
#************************************************************************************************

    ## Results and Responses:

        ### Abs response x_eq(t):
fig1 = plt.figure()
rect = fig1.patch

graphEQ = fig1.add_subplot(1,1,1)
graphEQ.plot(time, X, linewidth=2)
#graphEQ.plot(time, V, linewidth=2)
graphEQ.plot(time, SeismicY, linewidth=1)
#graphEQ.plot(time, SeismicV, linewidth=1)
graphEQ.set_title(' "SAW Seismic System Response" ', fontsize=25)
graphEQ.set_xlabel('Time [s]', fontsize=20)
graphEQ.legend(['x(t) [m]', 'y(t) [m]'], fontsize=13, loc='lower right')
graphEQ.grid(True)



plt.show()
