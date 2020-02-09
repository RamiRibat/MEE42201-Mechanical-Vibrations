# Mech Vib Tsak 2 (2 DOF, Wind Horizontal Vibration, Analytic Solution) #

#************************************************************************************************

    ## The required Libraries to be imported:

import numpy as np
from numpy import rad2deg
from numpy.linalg import inv, det
import matplotlib.pyplot as plt
from math import pi, sqrt, sin, cos, atan, exp
from scipy.linalg import eigh
from cmath import phase
import sympy as sp
from sympy import symbols
from sympy.physics.vector import init_vprinting, vlatex
init_vprinting(use_latex='mathjax', pretty_print=False)

#************************************************************************************************

    ## Functions:

def Z(w_w): # Impedences
    Z = np.matrix(np.complex_(np.zeros((4,4))))
    for i in range(0,4):
        for j in range(0,4):
            Z[i,j] = ( K[i,j]-w_w**2*M[i,j] + complex(0,w_w*C[i,j]) ) #+ complex(0,w_w*C[i,j]
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
m1 = rho * V # Floor mass [kg]
m2 = 1.5*m1

I = (l * w**4) / 12 # Column inertia [m^4]
ko = (3 * E * I) / H**3 # Column stiffness (cantilever beam)
k1 = 4*ko
k2 = k1

zeta = 0.05 # Daming ratio of the building
zeta_list = [0.05, 0.1, 0.15, 0.25, 0.5, 1, 2]
ze1, ze2, ze3, ze4, ze5, ze6, ze7 = symbols('zeta1, zeta2, zeta3, zeta4, zeta5, zeta6, zeta7')
ze_list = [ze1, ze2, ze3, ze4, ze5, ze6, ze7]
R = np.linspace(0.00001,4,100) # Frequency ratio list
w_list = np.linspace(0.00001,140.00001,140)
wn = sqrt(k1/m1) # Natural frequency
wd = sqrt(1-zeta**2)*wn # Damped Natural frequency

cc =  2*sqrt(m1*k1) # Critical damping
c1 = zeta * cc # Damping coefficient
c2 = c1

#************************************************************************************************

    ## Wind Input Data (fi = Fwio * cos(w_wi*t), direct force):

rho_a = 1.2 # Air density [kg/m^3]
A_w = L*H # Wind faced area of the building
v_w = 140 # Wind velocity [m/s]

#*# Excitation frequency [Hz]
freq_w1 = 20
freq_w2 = 20

#*# Cicular frequency [rad/s]
w_w1 = 2 * pi * freq_w1
w_w2 = 2 * pi * freq_w2
w_w = np.array((w_w1,w_w2))

#*# Wind force [N]
Fw1o = 0.5 * rho_a * v_w**2 * A_w 
Fw2o = 0.5 * rho_a * v_w**2 * A_w
Fwo = np.array((Fw1o,0,Fw2o,0))

#************************************************************************************************

     ## Absorber Data:

ma = 200
ka = ma*w_w1**2
ca = 0
mb = 200
kb = mb*w_w2**2
cb = 0
wna = sqrt(ka/ma)

M = np.array([[m1,0,0,0],[0,ma,0,0],[0,0,m2,0],[0,0,0,mb]])
C = np.array([[c1+c2+ca,-ca,-c2,0],[-ca,ca,0,0],[-c2,0,c2+cb,-cb],[0,0,-cb,cb]])
K = np.array([[k1+k2+ka,-ka,-k2,0],[-ka,ka,0,0],[-k2,0,k2+kb,-kb],[0,0,-kb,kb]])

        ### Natural modes & Natural frequencies:
evals, evecs = eigh(K, M)
wn1 = np.sqrt(evals[0])
wn2 = np.sqrt(evals[1])
wn3 = np.sqrt(evals[2])
wn4 = np.sqrt(evals[3])

#************************************************************************************************

    ## Simulation:

dt = 1/(8*freq_w1) # Time step
end_time = 1 # Finish time
time = np.arange(0, end_time, dt) # time range

        ### Initial Conditions:
x1o = 0.0
v1o = 0.0
x2o = 0.0
v2o = 0.0

        ### Response Histories & Amplitudes:
Fw1_history = [] # Wind input history
Fw2_history = [] # Wind input history
#-----------------------------------------------------------
Xw1_history = [] # Wind 1st response history
N1 = Z(w_w1)
for i in range(0,4):
    N1[i,0] = Fwo[i]
Xw1 = det(N1) / det(Z(w_w1))

Xwa_history = [] # Absorber response history
Na = Z(w_w1)
for i in range(0,4):
    Na[i,1] = Fwo[i]
Xwa = det(Na) / det(Z(w_w1))


Xw2_history = [] # Wind 2nd response history
N2 = Z(w_w2)
for i in range(0,4):
    N2[i,2] = Fwo[i]
Xw2 = det(N2) / det(Z(w_w2))

Xwb_history = [] # Absorber response history
Nb = Z(w_w2)
for i in range(0,4):
    Nb[i,3] = Fwo[i]
Xwb = det(Nb) / det(Z(w_w2))

        ### Time Responses:
for t in time:
#--------------------------------------------------
    fw1 = Fw1o * cos(w_w1*t) # input wind force
    Fw1_history.append(0.00000001*fw1)
    
    fw2 = Fw2o * cos(w_w2*t) # input wind force
    Fw2_history.append(0.00000001*fw2)
#--------------------------------------------------
    #* 1st steady state response:    
    x_w1 = abs(Xw1)*cos(w_w1*t-phase(Xw1)) # Steady state response from wind
    Xw1_history.append(x_w1*1e2)

    #* Absorber a steady state response:
    x_wa = abs(Xwa)*cos(w_w1*t-phase(Xwa)) # Steady state response from wind
    Xwa_history.append(x_wa*1e2)

    #* 2nd steady state response:
    x_w2 = abs(Xw2)*cos(w_w2*t-phase(Xw2)) # Steady state response from wind        
    Xw2_history.append(x_w2*1e2)
    
    #* Absorber b steady state response:
    x_wb = abs(Xwb)*cos(w_w2*t-phase(Xwb)) # Steady state response from wind
    Xwb_history.append(x_wb*1e2)
#--------------------------------------------------
Xw1_w_history = np.zeros((len(w_list),1)) # Earhquacke absolute response history
Xwa_w_history = np.zeros((len(w_list),1)) # Earhquacke absolute response history
Xw2_w_history = np.zeros((len(w_list),1)) # Earhquacke absolute response history
Xwb_w_history = np.zeros((len(w_list),1)) # Earhquacke absolute response history
for j in range(0,len(w_list)):
        w = w_list[j]
        
        N1 = Z(w)
        for i in range(0,4):
            N1[i,0] = Fwo[i]
        Xw1 = det(N1) / det(Z(w))

        Na = Z(w)
        for i in range(0,4):
            Na[i,1] = Fwo[i]
        Xwa = det(Na) / det(Z(w))
        
        N2 = Z(w)
        for i in range(0,4):
            N2[i,2] = Fwo[i]
        Xw2 = det(N2)/det(Z(w))

        Nb = Z(w)
        for i in range(0,4):
            Nb[i,3] = Fwo[i]
        Xwb = det(Nb)/det(Z(w))
        
        Xw1_w_history[j] = abs(Xw1*1e2)
        Xwa_w_history[j] = abs(Xwa*1e2)
        Xw2_w_history[j] = abs(Xw2*1e2)
        Xwb_w_history[j] = abs(Xwb*1e2)

#************************************************************************************************

    ## Results and Responses:

        ### Absolute Response
fig1 = plt.figure()
rect = fig1.patch

graphW = fig1.add_subplot(1,1,1)
graphW.plot(time, Xw1_history,  linewidth=2)
graphW.plot(time, Xwa_history,  linewidth=2)
graphW.plot(time, Xw2_history,  linewidth=2)
graphW.plot(time, Xwb_history,  linewidth=2)
#graphW.plot(time, Fw1_history,  linewidth=1.5)
#graphW.plot(time, Fw2_history,  linewidth=1.5)
graphW.set_title('"System Responses for Wind"', fontsize=25)
graphW.set_xlabel('Time [s]', fontsize=20)
graphW.legend(['x_w1(t) [cm]', 'x_wa(t) [cm]', 'x_w2(t) [cm]', 'x_wb(t) [cm]'], fontsize=13, loc='lower right')
#graphW.legend(['f_w1(t) [m]', 'f_w2(t) [m]'], fontsize=10, loc='lower right')
graphW.grid(True)


        ### Frequency Response
fig2 = plt.figure()
rect = fig2.patch

graphW2 = fig2.add_subplot(1,1,1)
graphW2.plot(w_list, Xw1_w_history, color='pink', linewidth=10)
graphW2.plot(w_list, Xw2_w_history, color='lightgreen', linewidth=5)
graphW2.plot(w_list, Xwa_w_history, color='darkred', linewidth=2)
graphW2.plot(w_list, Xwb_w_history, color='darkgreen', linewidth=2)
graphW2.set_title('"Floors Amplitudes with w_w"', fontsize=25)
graphW2.set_ylabel('X', fontsize=20)
graphW2.set_xlabel('w_w [rad/s]', fontsize=20)
graphW2.legend(['Xw1(t) [cm]', 'Xw2(t) [cm]', 'Xwa(t) [cm]', 'Xwb(t) [cm]'], fontsize=13, loc='upper right')
graphW2.grid(True)


plt.show()
