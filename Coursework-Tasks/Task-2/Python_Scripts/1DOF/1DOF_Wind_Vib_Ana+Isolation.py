# Mech Vib Tsak 2 (SDOF, Wind Vibration, Analytic Sol) #

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
zeta = 10 # Daming ratio of the building
zeta_list = [0.05, 0.1, 0.15, 0.25, 0.5, 1, 2]
ze1, ze2, ze3, ze4, ze5, ze6, ze7 = symbols('zeta1, zeta2, zeta3, zeta4, zeta5, zeta6, zeta7')
ze_list = [ze1, ze2, ze3, ze4, ze5, ze6, ze7]
R = np.linspace(0.00001,4.00001,120) # Frequency ratio list
wn = sqrt(k/m) # Natural frequency
#wd = sqrt(1-zeta**2)*wn # Damped Natural frequency
cc = 2*sqrt(m*k) # Critical damping
c = zeta * cc # Damping coefficient

#************************************************************************************************

    ## Input Data (Sin wave):

        ### Wind (f = Fo * cos(w_w*t), direct force):
freq_w = 20 # Excitation frequency [Hz]
w_w = 2 * pi * freq_w # Cicular frequency [rad/s]
v_w = 140 # Wind velocity [m/s]
rho_a = 1.2 # Air density [kg/m^3]
A_w = L*H # Wind faced area of the building 
Fwo = 0.5 * rho_a * v_w**2 * A_w # Wind force [N]

#************************************************************************************************

    ## Simulation:

        ### Time Installation:
dt = 1/(16*freq_w) # Time step
end_time = 0.5 # Finish time
time = np.arange(0, end_time, dt) # Time range

        ### Initial Conditions:
xo = 0.0
vo = 0.0

        ### Response Histories & Amplitudes:
r_w = w_w/wn

Fw_history = [] # Wind input history

x_st = Fwo/k # Static defliction

Xw_history = np.zeros((len(zeta_list),len(time))) # Wind response history

        ### Time Responses:
#for t in time:
#--------------------------------------------------
    #fw = Fwo * cos(w_w*t) # input wind force
    #Fw_history.append(fw*1e-6) # MN
#--------------------------------------------------
for i in range(0,len(zeta_list)):
    ze = zeta_list[i]
    c = ze * cc # Damping coefficient
    Xw = Fwo / ((k-m*w_w**2) + complex(0,w_w*c))
    for j in range(0,len(time)):
        t = time[j]
        x_w = abs(Xw)*cos(w_w*t-phase(Xw)) # Steady state response from wind
        Xw_history[i,j] = x_w*1e6 # micro
    
#************************************************************************************************

    ## Results and Responses:
fig1 = plt.figure()
rect = fig1.patch

for i in range(0,len(zeta_list)):
    graphW = fig1.add_subplot(1,1,1)
    graphW.plot(time, Xw_history[i,0:], linewidth=2.0)
    #graphW.plot(time, Fw_history, linewidth=1.0)
    #graphW.legend(['x_w(t) [cm]', 'f(t) [MN]'], fontsize=12, loc='lower right')
    graphW.set_title(' "Isolation for Wind at resonace, [microm]" ', fontsize=25)
    graphW.legend( ["${}$".format( vlatex(ze_list[k])+' = '+str(zeta_list[k])) for k in range(0,len(zeta_list))  ] , fontsize=13, loc='upper right' )
    graphW.set_xlabel('Time [s]', fontsize=20)
    graphW.grid(True)

plt.show()
