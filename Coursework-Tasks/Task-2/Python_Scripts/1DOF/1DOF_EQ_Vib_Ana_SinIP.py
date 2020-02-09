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
zeta = 0.05 # Daming ratio of the building
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
freq_eq = 20 # Excitation frequency [Hz]
w_eq = 2 * pi * freq_eq # Cicular frequency [rad/s]
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

Xeq_history = [] # Earhquacke absolute response history
Xeq = Yo * (k+complex(0,w_eq*c)) / ((k-m*w_eq**2)+complex(0,w_eq*c))

Q_history = [] # Earhquacke relative response history
Q = (Yo*m*w_eq**2) / (k-m*w_eq**2+complex(0,c))

        ### Time Responses:
for t in time:
#--------------------------------------------------
    y = Yo * sin(w_eq*t) # input seismic force
    Y_history.append(y)
#--------------------------------------------------
    x_eq = abs(Xeq)*sin(w_eq*t-phase(Xeq)) ## Steady state response from earthquack
    Xeq_history.append(x_eq)
#--------------------------------------------------
Td_history = np.zeros((len(zeta_list),len(R))) # Displacement transmissibilty
TH_history = np.zeros((len(zeta_list),len(R))) # Phase
for i in range(0,len(zeta_list)):
    for j in range(0,len(R)):
        ze = zeta_list[i]
        r = R[j]
        Td = sqrt( 1 + (2*ze*r)**2 ) / sqrt( (1-r**2)**2 + (2*ze*r)**2 )
        Td_history[i,j] = Td
        th = rad2deg(atan(2*ze*r**3 / (1+(4*ze**2-1)*r**2)))
        if th >= 0:
            TH = rad2deg(atan(2*ze*r**3 / (1+(4*ze**2-1)*r**2)))
        elif th < 0:
            TH = 180+rad2deg(atan(2*ze*r**3 / (1+(4*ze**2-1)*r**2)))
        else:
            TH = 0
        TH_history[i,j] = TH
#--------------------------------------------------
Ft_kY_history = np.zeros((len(zeta_list),len(R))) # Force transmitted ratio
for i in range(0,len(zeta_list)):
    for j in range(0,len(R)):
        ze = zeta_list[i]
        r = R[j]
        Tf = ( r**2 * sqrt( 1 + (2*ze*r)**2 ) ) / sqrt( (1-r**2)**2 + (2*ze*r)**2 )
        Ft_kY_history[i,j] = Tf
#-------------------------------------------------
for t in time:
    q = abs(Q)*sin(w_eq*t-phase(Q))
    Q_history.append(q)
    
#************************************************************************************************

    ## Results and Responses:

        ### Abs response x_eq(t):
fig1 = plt.figure()
rect = fig1.patch

graphEQ = fig1.add_subplot(1,1,1)
graphEQ.plot(time, Xeq_history, linewidth=2.0)
graphEQ.plot(time, Y_history, linewidth=1.0)
graphEQ.legend(['x_eq(t) [m]', 'y(t) [m]'], fontsize=13, loc='lower right')
graphEQ.set_title('Seismic Absolute Response', fontsize=25)
graphEQ.set_xlabel('Time [s]', fontsize=20)
graphEQ.grid(True)


        ### Displacement transmissibilty
fig2a = plt.figure()
rect = fig2a.patch
fig2b = plt.figure()
rect = fig2b.patch
for i in range(0,len(zeta_list)):
    graphEQ2a = fig2a.add_subplot(1,1,1)
    graphEQ2a.plot(R, Td_history[i,0:], linewidth=1.5)
    graphEQ2a.legend( ["${}$".format( vlatex(ze_list[k])+' = '+str(zeta_list[k])) for k in range(0,len(zeta_list))  ] , fontsize=13, loc='upper right' )
    graphEQ2a.set_title('Displacement transmissibilty with r for Earthquake', fontsize=25)
    graphEQ2a.set_ylabel('Td = X/Y', fontsize=20)
    graphEQ2a.set_xlabel('r = w_eq/wn', fontsize=20)
    graphEQ2a.grid(True)
    graphEQ2b = fig2b.add_subplot(1,1,1)
    graphEQ2b.plot(R, TH_history[i,0:], linewidth=1.5)
    graphEQ2b.legend( ["${}$".format( vlatex(ze_list[k])+' = '+str(zeta_list[k])) for k in range(0,len(zeta_list))  ] , fontsize=13, loc='lower right' )
    graphEQ2b.set_title('Phase angle with r for Earthquake', fontsize=25)
    graphEQ2b.set_ylabel('phase [deg]', fontsize=20)
    graphEQ2b.set_xlabel('r = w_eq/wn', fontsize=20)
    graphEQ2b.grid(True)


        ### Force transmitted:
fig3 = plt.figure()
rect = fig3.patch
for i in range(0,len(zeta_list)):
    graphEQ3 = fig3.add_subplot(1,1,1)
    graphEQ3.plot(R, Ft_kY_history[i,0:], linewidth=1.5)
    graphEQ3.legend( ["${}$".format( vlatex(ze_list[k])+' = '+str(zeta_list[k])) for k in range(0,len(zeta_list))  ] , fontsize=13, loc='upper left' )
    graphEQ3.set_title('Force transmitted with r for Earthquake', fontsize=25)
    graphEQ3.set_ylabel('Tf = Ft/kYo', fontsize=20)
    graphEQ3.set_xlabel('r = w_eq/wn', fontsize=20)
    graphEQ3.grid(True)


        ### Relative response q(t)
fig4 = plt.figure()
rect = fig4.patch

graphEQ4 = fig4.add_subplot(1,1,1)
graphEQ4.plot(time, Q_history, linewidth=2.0)
graphEQ4.plot(time, Y_history, linewidth=1.0)
graphEQ4.set_title('Relative response for Earthquake', fontsize=25)
graphEQ4.legend(['q(t) [m]', 'y(t) [m]'], fontsize=13, loc='lower right')
graphEQ4.set_xlabel('Time [s]', fontsize=20)
graphEQ4.grid(True)


plt.show()
