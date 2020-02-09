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
m1 = rho * V # Floor mass [kg]
m2 = m1

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

M = np.matrix([[m1,0],
               [0,m2]])
C = np.array([[c1+c2,-c2],
              [-c2,c2]])
K = np.array([[k1+k2,-k2],
              [-k2,k2]])

                #### Natural modes & Natural frequencies:
evals, evecs = eigh(K, M)
wn1 = np.sqrt(evals[0])
wn2 = np.sqrt(evals[1])

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
Fwo = sp.Matrix([Fw1o, Fw2o])
#-----------------------------------------------------------------

    ## Earth quicks Input Data (y = Yo * sin(w_eq*t, base motion)):
freq_eq = 20 # Excitation frequency [Hz]
w_eq = 2 * pi * freq_eq # Cicular frequency [rad/s]
Yo = 5 # Earth quick amplitude [m]
F1o = w_eq**2*m1*Yo
F2o = w_eq**2*m2*Yo
Fo = sp.Matrix([F1o, F2o]) # Equivalent Seismic Forces

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
#--------------------------------------------------
Fw1_history = [] # Wind input history
Fw2_history = [] # Wind input history

Xw1_history = [] # Wind 1st response history
N1 = Z(w_w1)
for i in range(0,2):
    N1[i,0] = Fwo[i]
Xw1 = det(N1)/det(Z(w_w1))

Xw2_history = [] # Wind 2nd response history
N2 = Z(w_w2)
for i in range(0,2):
    N2[i,1] = Fwo[i]
Xw2 = det(N2)/det(Z(w_w2))
#--------------------------------------------------
Y_history = [] # Earh quick input history

Xeq1_history = [] # Earhquack m1 response history
N1 = Z(w_eq)
N1[0,0] = (complex(k1, w_eq*c1))*Yo
Xeq1 = det(N1)/det(Z(w_eq))

Xeq2_history = [] # Earhquack m2 response history
N2 = Z(w_eq)
N2[1,0] = (complex(k1, w_eq*c1))*Yo
Xeq2 = det(N2)/det(Z(w_eq))
#--------------------------------------------------
#--------------------------------------------------
for t in time:
#------------------------------------------------------------------------------------------
    fw1 = Fw1o * cos(w_w1*t) # input wind force
    Fw1_history.append(0.00000001*fw1)
    
    fw2 = Fw2o * cos(w_w2*t) # input wind force
    Fw2_history.append(0.00000001*fw2)


    y = Yo * sin(w_eq*t)
    Y_history.append(y)
#------------------------------------------------------------------------------------------
    #* 1st steady state response:    
    x_w1 = abs(Xw1)*cos(w_w1*t-phase(Xw1))
    Xw1_history.append(x_w1*1e2)

    #* 2nd state response
    x_w2 = abs(Xw2)*cos(w_w2*t-phase(Xw2))
    Xw2_history.append(x_w2*1e2)


            #* 1st steady state response:    
    xeq1 = abs(Xeq1)*sin(w_eq*t-phase(Xeq1))  # Steady state response from wind
    Xeq1_history.append(xeq1)
    #* 2nd state response
    xeq2 = abs(Xeq2)*sin(w_eq*t-phase(Xeq2)) # Steady state response from wind
    Xeq2_history.append(xeq2)
#------------------------------------------------------------------------------------------
T1_history = np.zeros((len(time),1))
T2_history = np.zeros((len(time),1))
for i in range(0,len(time)):
    T1_history[i] = Xw1_history[i] + Xeq1_history[i] # Total Steady state response
    T2_history[i] = Xw2_history[i] + Xeq2_history[i] # Total Steady state response

#************************************************************************************************

    ## Results and Responses:

        ### Absolute Response:
fig1 = plt.figure()
rect = fig1.patch

graphT = fig1.add_subplot(1,1,1)
graphT.plot(time, Xw1_history, linewidth=2.0)
graphT.plot(time, Xw2_history, linewidth=2.0)
graphT.plot(time, Fw1_history, linewidth=1.0)
graphT.plot(time, Fw2_history, linewidth=1.0)
graphT.plot(time, Xeq1_history, linewidth=2.0)
graphT.plot(time, Xeq2_history, linewidth=2.0)
graphT.plot(time, Y_history, linewidth=1.0)
graphT.plot(time, T1_history, linewidth=2.0)
graphT.plot(time, T2_history, linewidth=2.0)
graphT.legend(['x_w1(t) [m]', 'x_w1(t) [m]', 'fw1(t) [MN]', 'fw2(t) [MN]', 'x_eq1(t) [m]', 'x_eq2(t) [m]', 'y(t) [m]', 'Total [m]'], fontsize=13, loc='lower right')
graphT.set_title(' "Total System Response" ', fontsize=25)
graphT.set_xlabel('Time [s]', fontsize=20)
graphT.grid(True)

plt.show()
