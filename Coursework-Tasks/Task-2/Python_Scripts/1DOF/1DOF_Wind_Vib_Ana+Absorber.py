# Mech Vib Tsak 2 (Single DOF, Wind Vibration + Absorber, Analytic Solution):

#************************************************************************************************

    ## The required Libraries to be imported:

import numpy as np
from numpy import rad2deg
from numpy.linalg import inv, det
import matplotlib.pyplot as plt
from math import pi, sqrt, sin, cos, atan
from cmath import phase
from sympy import symbols
from scipy.linalg import eigh
from sympy.physics.vector import init_vprinting, vlatex
init_vprinting(use_latex='mathjax', pretty_print=False)


#************************************************************************************************

    ## Functions:

def Z(w): # Impedences
    Z = np.matrix(np.complex_(np.zeros((2,2))))
    for i in range(0,2):
        for j in range(0,2):
            Z[i,j] = ( K[i,j]-w**2*M[i,j] + complex(0,w*C[i,j]) ) #+ complex(0,w_w*C[i,j]
    return Z

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
k = ko # Equivelant columns stiffnesses
zeta = 0.05 # Daming ratio of the building
zeta_list = [0.05, 0.1, 0.15, 0.25, 0.5, 1, 2]
ze1, ze2, ze3, ze4, ze5, ze6, ze7 = symbols('zeta1, zeta2, zeta3, zeta4, zeta5, zeta6, zeta7')
ze_list = [ze1, ze2, ze3, ze4, ze5, ze6, ze7]
R = np.linspace(0.00001,4,100) # Frequency ratio list
w_list = np.linspace(0.00001,140.00001,140)
wn = sqrt(k/m) # Natural frequency
wd = sqrt(1-zeta**2)*wn # Damped Natural frequency
cc = 2*sqrt(m*k) # Critical damping
c = zeta * cc # Equivelant damping coefficient

#************************************************************************************************

    ## Input Data (Sin wave):

        ### Wind (f = Fo * cos(w_w*t), direct force):
freq_w = 20 # Excitation frequency [Hz]
w_w =  2 * pi * freq_w # Cicular frequency [rad/s]
v_w = 140 # Wind velocity [m/s]
rho_a = 1.2 # Air density [kg/m^3]
A_w = L*H # Wind faced area of the building 
Fwo = 0.5 * rho_a * v_w**2 * A_w # Wind force [N]


#************************************************************************************************

     ## Absorber Data:

ma = 200
ka = ma*w_w**2
ca = 0
#wna = sqrt(ka/ma)

M = np.matrix([[m,0], [0,ma]])
C = np.array([[c+ca,-ca], [-ca,ca]])
K = np.array([[k+ka,-ka], [-ka,ka]])

evals, evecs = eigh(K,M)
wn1 = np.sqrt(evals[0])
wn2 = np.sqrt(evals[1])

#************************************************************************************************


    ## Simulation:

dt = 1/(8*freq_w) # Time step
end_time = 0.5 # Finish time
time = np.arange(0, end_time, dt) # time range

        ### Initial Conditions:
xo = 0.0
vo = 0.0

        ### Response Histories & Amplitudes:
r_w = w_w/wn
#r_wa = w_w/wna

Fw_history = [] # Wind input history

Xw_history = [] # Wind response history
Xw = (Z(w_w)[1,1]*Fwo) / ( (Z(w_w)[0,0]*Z(w_w)[1,1]) - (Z(w_w)[0,1])**2 )

Xwa_history = [] # Absorber response history
Xwa = - (Z(w_w)[0,1]*Fwo) / ( (Z(w_w)[0,0]*Z(w_w)[1,1]) - (Z(w_w)[0,1])**2 )


        ### Time Responses:
for t in time:
#--------------------------------------------------
    fw = Fwo * cos(w_w*t) # input wind force
    Fw_history.append(fw*1e-3)
#--------------------------------------------------
    x_w = abs(Xw)*cos(w_w*t-phase(Xw))
    Xw_history.append(x_w*1e2)
    
    x_wa = abs(Xwa)*cos(w_w*t-phase(Xwa))
    Xwa_history.append(x_wa*1e2)
#------------------------------------------------------------------------------------------
Xw_w_history = np.zeros((len(w_list),1)) # Earhquacke absolute response history
Xwa_w_history = np.zeros((len(w_list),1)) # Earhquacke absolute response history
for j in range(0,len(w_list)):
        w = w_list[j]
        Xw = (Z(w)[1,1]*Fwo) / ( (Z(w)[0,0]*Z(w)[1,1]) - (Z(w)[0,1])**2 )
        Xwa = - (Z(w)[0,1]*Fwo) / ( (Z(w)[0,0]*Z(w)[1,1]) - (Z(w)[0,1])**2 )
        Xw_w_history[j] = abs(Xw)
        Xwa_w_history[j] = abs(Xwa)


#************************************************************************************************

    ## Results and Responses:
fig1 = plt.figure()
rect = fig1.patch

graphW = fig1.add_subplot(1,1,1)
graphW.plot(time, Xw_history, linewidth=2.0)
graphW.plot(time, Xwa_history, linewidth=2.0)
graphW.plot(time, Fw_history, linewidth=1.0)
graphW.legend(['x_w(t) [cm]', 'x_wa(t) [cm]', 'f(t) [kN]'], fontsize=12, loc='lower right')
graphW.set_title(' "System Response" ', fontsize=25)
graphW.set_xlabel('Time [s]', fontsize=20)
graphW.grid(True)

fig2 = plt.figure()
rect = fig2.patch
graphw2 = fig2.add_subplot(1,1,1)
graphw2.plot(w_list, Xw_w_history, color='orange', linewidth=6)
graphw2.plot(w_list, Xwa_w_history, color='darkblue', linewidth=2)
graphw2.set_title('"Floor & Absorber Amplitude"', fontsize=25)
graphw2.set_ylabel('X ', fontsize=20)
graphw2.set_xlabel('w_w', fontsize=20)
graphw2.legend(['Xw(t) [m]', 'Xwa(t) [m]'], fontsize=13, loc='upper right')
graphw2.grid(True)


plt.show()
