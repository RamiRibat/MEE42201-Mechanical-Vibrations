# Mech Vib Tsak 2 (SDOF, Earthquake Vibration, Analytic Sol) #

#************************************************************************************************

    ## The required Libraries to be imported:

import numpy as np
from numpy import rad2deg
from numpy.linalg import inv, det
import matplotlib.pyplot as plt
from math import pi, sqrt, sin, cos, atan
from cmath import phase
from sympy import symbols
from sympy.physics.vector import init_vprinting, vlatex
init_vprinting(use_latex='mathjax', pretty_print=False)


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
zeta = 0.9 # Daming ratio of the building
zeta_list = [0.05, 0.1, 0.15, 0.25, 0.5, 1, 2]
ze1, ze2, ze3, ze4, ze5, ze6, ze7 = symbols('zeta1, zeta2, zeta3, zeta4, zeta5, zeta6, zeta7')
ze_list = [ze1, ze2, ze3, ze4, ze5, ze6, ze7]
R = np.linspace(0.00001,4,100) # Frequency ratio list
wn = sqrt(k/m) # Natural frequency
wd = sqrt(1-zeta**2)*wn # Damped Natural frequency
cc = 2*sqrt(m*k) # Critical damping
c = zeta * cc # Equivelant damping coefficient


#************************************************************************************************

    ## Input Data (Sin wave):

        #### Earth quicks (y = Yo * sin(w_eq*t, base motion)):
freq_eq = 18 # Excitation frequency [Hz]
w_eq = 64# 2 * pi * freq_eq # Cicular frequency [rad/s]
Yo = 5 # Earth quick amplitude [m]


#************************************************************************************************

    ## Simulation:

dt = 1/(8*freq_eq) # Time step
end_time = 1 # Finish time
time = np.arange(0, end_time, dt) # time range

        ### Initial Conditions:
xo = 0.0
vo = 0.0

        ### Response Histories & Amplitudes:
Y_history = [] # Earh quick input history
r_eq = w_eq/wn

Xeq_history = np.zeros((len(zeta_list),len(time))) # Earhquacke absolute response history

        ### Time Responses:
#for t in time:
#--------------------------------------------------
    #y = Yo * sin(w_eq*t)
    #Y_history.append(y)
#--------------------------------------------------
for i in range(0,len(zeta_list)):
    ze = zeta_list[i]
    c = ze * cc # Damping coefficient
    Xeq = Yo * (k+complex(0,w_eq*c)) / ((k-m*w_eq**2)+complex(0,w_eq*c))
    for j in range(0,len(time)):
        t = time[j]
        x_eq = abs(Xeq)*sin(w_eq*t-phase(Xeq)) #-phase(Xeq)
        Xeq_history[i,j] = x_eq*1e2 # cm
    
#************************************************************************************************

    ## Results and Responses:

        ### Abs response x_eq(t):
fig1 = plt.figure()
rect = fig1.patch

for i in range(0,len(zeta_list)):
    graphEQ = fig1.add_subplot(1,1,1)
    graphEQ.plot(time, Xeq_history[i,0:], linewidth=2.0)
    #graphEQ.plot(time, Y_history, linewidth=1.0)
    graphEQ.legend( ["${}$".format( vlatex(ze_list[k])+' = '+str(zeta_list[k])) for k in range(0,len(zeta_list))  ] , fontsize=13, loc='upper right' )
    graphEQ.set_title(' "Isolation for EQ at resonance,[cm]"', fontsize=25)
    graphEQ.set_xlabel('Time [s]', fontsize=20)
    graphEQ.grid(True)




plt.show()
