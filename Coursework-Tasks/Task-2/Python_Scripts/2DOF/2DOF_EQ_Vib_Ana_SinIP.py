# Mech Vib Tsak 2 (2DOF, WEarthquake Vibration, Analytic Sol) #

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

    ## Earth quicks Input Data (y = Yo * sin(w_eq*t, base motion)):
freq_eq = 20 # Excitation frequency [Hz]
w_eq = 2 * pi * freq_eq # Cicular frequency [rad/s]
Yo = 5 # Earth quick amplitude [m]
F1o = w_eq**2*m1*Yo
F2o = w_eq**2*m2*Yo
Fo = sp.Matrix([F1o, F2o]) # Equivalent Seismic Forces


#************************************************************************************************

    ## Simulation:

dt = 1/(8*freq_eq) # Time step
end_time = 1 # Finish time
time = np.arange(0, end_time, dt) # time range

        ### Initial Conditions:
x1o = 0.0
v1o = 0.0
x2o = 0.0
v2o = 0.0

        ### Response Histories & Amplitudes:
Y_history = [] # Earh quick input history
#-----------------------------------------------------------
Xeq1_history = [] # Earhquack m1 response history
N1 = Z(w_eq)
N1[0,0] = (complex(k1, w_eq*c1))*Yo
Xeq1 = det(N1)/det(Z(w_eq))

Xeq2_history = [] # Earhquack m2 response history
N2 = Z(w_eq)
N2[1,0] = (complex(k1, w_eq*c1))*Yo
Xeq2 = det(N2)/det(Z(w_eq))
#-----------------------------------------------------------
Q1_history = [] # Earhquack m1 relative response history
Nq1 = Z(w_eq)
for i in range(0,2):
    Nq1[i,0] = Fo[i]
Q1 = det(N1)/det(Z(w_eq))

Q2_history = [] # Earhquack m2 relative response history
Nq2 = Z(w_eq)
for i in range(0,2):
    Nq2[i,1] = Fo[i]
Q2 = det(N2)/det(Z(w_eq))

        ### Time Response:
for t in time:
#--------------------------------------------------
    y = Yo * sin(w_eq*t)
    Y_history.append(y)
#--------------------------------------------------
                #### Absolute response x_eq = xeq1 + xeq2
    #* 1st steady state response:    
    xeq1 = abs(Xeq1)*sin(w_eq*t-phase(Xeq1)) # +  abs(Xeq12)*cos(w_eq*t-phase(Xeq12)) # Steady state response from wind
    Xeq1_history.append(xeq1)
    #* 2nd state response
    xeq2 = abs(Xeq2)*sin(w_eq*t-phase(Xeq2)) # + abs(Xeq22)*cos(w_eq*t-phase(Xeq22))# Steady state response from wind
    Xeq2_history.append(xeq2)
#--------------------------------------------------
            #### Relative response q = x_eq - y 
    #* 1st steady state response:    
    q1 = abs(Q1)*sin(w_eq*t-phase(Q1))  # Steady state response from wind
    Q1_history.append(q1)
    #* 2nd state response
    q2 = abs(Q2)*sin(w_eq*t-phase(Q2))  # Steady state response from wind
    Q2_history.append(q2)
#--------------------------------------------------
Xeq1_w_history = np.zeros((len(w_list),1)) # Earhquacke absolute response history
Xeq2_w_history = np.zeros((len(w_list),1)) # Earhquacke absolute response history
for j in range(0,len(w_list)):
        w = w_list[j]
        
        N1 = Z(w)
        for i in range(0,2):
            N1[i,0] = Fo[i]
        Xeq1 = det(N1)/det(Z(w))
        
        N2 = Z(w)
        for i in range(0,2):
            N2[i,1] = Fo[i]
        Xeq2 = det(N2)/det(Z(w))
        
        Xeq1_w_history[j] = abs(Xeq1*1e2)
        Xeq2_w_history[j] = abs(Xeq2*1e2)

#************************************************************************************************

    ## Results and Responses:

        ### Absolute Response:
fig1 = plt.figure()
rect = fig1.patch

graphEQ = fig1.add_subplot(1,1,1)
graphEQ.plot(time, Xeq1_history,  linewidth=2)
graphEQ.plot(time, Xeq2_history,  linewidth=2)
graphEQ.plot(time, Y_history,  linewidth=1)
graphEQ.set_title('"Absolute Response for Seismic"', fontsize=25)
graphEQ.set_xlabel('Time [s]', fontsize=20)
graphEQ.legend(['x_eq1(t) [m]', 'x_eq2(t) [m]', 'y(t) [m]'], fontsize=10, loc='lower right')
graphEQ.grid(True)


        ### Relative Response:
fig2 = plt.figure()
rect = fig2.patch

graphEQ2 = fig2.add_subplot(1,1,1)
graphEQ2.plot(time, Q1_history,  linewidth=2)
graphEQ2.plot(time, Q2_history,  linewidth=2)
graphEQ2.plot(time, Y_history,  linewidth=1)
graphEQ2.set_title('"Relative Response"', fontsize=25)
graphEQ2.set_xlabel('Time [s]', fontsize=20)
graphEQ2.legend(['q1(t) [m]', 'q2(t) [m]', 'y(t) [m]'], fontsize=10, loc='lower right')
graphEQ2.grid(True)

        ### Frequency Response
fig3 = plt.figure()
rect = fig3.patch
graphEQ3 = fig3.add_subplot(1,1,1)
graphEQ3.plot(w_list, Xeq1_w_history, linewidth=2)
graphEQ3.plot(w_list, Xeq2_w_history, linewidth=2)
graphEQ3.set_title('"Floors Amplitudes with w_w"', fontsize=25)
graphEQ3.set_ylabel('X', fontsize=20)
graphEQ3.set_xlabel('w_w [rad/s]', fontsize=20)
graphEQ3.legend(['Xw1(t) [cm]', 'Xw2(t) [cm]'], fontsize=13, loc='upper right')
graphEQ3.grid(True)

plt.show()
