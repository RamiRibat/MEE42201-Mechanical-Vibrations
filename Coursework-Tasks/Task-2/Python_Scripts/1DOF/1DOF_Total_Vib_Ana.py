# Mech Vib Tsak 2 (SDOF, Total Vibration, Analytic Sol) #

#************************************************************************************************

    ## The required Libraries to be imported:

import numpy as np
from numpy import rad2deg
from numpy.linalg import inv
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
R = np.linspace(0.00001,4.00001,120) # Frequency ratio list
wn = sqrt(k/m) # Natural frequency
wd = sqrt(1-zeta**2)*wn # Damped Natural frequency
cc = 2*sqrt(m*k) # Critical damping
c = zeta * cc # Damping coefficient

#************************************************************************************************

    ## Input Data (Sin wave):

        ### Wind (f = Fo * cos(w_w*t), direct force):
freq_w = 10 # Excitation frequency [Hz]
w_w = 2 * pi * freq_w # Cicular frequency [rad/s]
v_w = 140 # Wind velocity [m/s]
rho_a = 1.2 # Air density [kg/m^3]
A_w = L*H # Wind faced area of the building 
Fwo = 0.5 * rho_a * v_w**2 * A_w # Wind force [N]

        ### Earthquacks (y = Yo * sin(w_eq*t), base motion)):
freq_eq = 20 # Excitation frequency [Hz]
w_eq = 2 * pi * freq_eq # Cicular frequency [rad/s]
Yo = 5 # Earth quick amplitude [m]

#************************************************************************************************

    ## Simulation:

        ### Time Installation:
dt = 1/(8*freq_w) # Time step
end_time = 1 # Finish time
time = np.arange(0, end_time, dt) # Time range

        ### Initial Conditions:
xo = 0.0
vo = 0.0

        ### Response Histories & Amplitudes:
#--------------------------------------------------
r_w = w_w/wn
Fw_history = [] # Wind input history
x_st = Fwo/k # Static defliction
Xw_history = [] # Wind response history
Xw = Fwo / ((k-m*w_w**2) + complex(0,w_w*c))
#--------------------------------------------------
r_eq = w_eq/wn
Y_history = [] # Earh quick input history
Xeq_history = [] # Earhquacke absolute response history
Xeq = Yo * (k+complex(0,w_eq*c)) / ((k-m*w_eq**2)+complex(0,w_eq*c))
#--------------------------------------------------
#--------------------------------------------------
        ### Time Responses:
for t in time:
#--------------------------------------------------
    fw = Fwo * cos(w_w*t) # input wind force
    Fw_history.append(fw*1e-6) # MN

    y = Yo * sin(w_eq*t) # input seismic force
    Y_history.append(y)
#--------------------------------------------------
    x_w = abs(Xw)*cos(w_w*t-phase(Xw)) # Steady state response from wind
    Xw_history.append(x_w) 

    x_eq = abs(Xeq)*sin(w_eq*t-phase(Xeq)) # Steady state response from earthquack
    Xeq_history.append(x_eq)
#--------------------------------------------------   
T_history = np.zeros((len(time),1))
for i in range(0,len(time)):
    T_history[i] = Xw_history[i] + Xeq_history[i] # Total Steady state response

#************************************************************************************************

    ## Results and Responses:

        ### Abs response x_eq(t):

fig1 = plt.figure()
rect = fig1.patch

graphT = fig1.add_subplot(1,1,1)
graphT.plot(time, Xw_history, linewidth=2.0)
graphT.plot(time, Fw_history, linewidth=1.0)
graphT.plot(time, Xeq_history, linewidth=2.0)
graphT.plot(time, Y_history, linewidth=1.0)
graphT.plot(time, T_history, linewidth=2.0)
graphT.legend(['x_w(t) [m]', 'f(t) [MN]', 'x_eq(t) [m]', 'y(t) [m]', 'Total [m]'], fontsize=12, loc='lower right')
graphT.set_title(' "Total System Response" ', fontsize=25)
graphT.set_xlabel('Time [s]', fontsize=20)
graphT.grid(True)

plt.show()
