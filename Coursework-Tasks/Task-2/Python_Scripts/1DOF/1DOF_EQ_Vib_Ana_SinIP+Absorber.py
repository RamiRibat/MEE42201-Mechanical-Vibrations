# Mech Vib Tsak 2 (SDOF, Earthquake Vibration + ABSORBER, Analytic Sol) #

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

        #### Earth quicks (y = Yo * sin(w_eq*t, base motion)):
freq_eq = 20 # Excitation frequency [Hz]
w_eq = 2 * pi * freq_eq # Cicular frequency [rad/s]
Yo = 5 # Earth quick amplitude [m]

#************************************************************************************************

     ## Absorber Design Data:

ma = 200
ka = ma*w_eq**2
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

dt = 1/(8*freq_eq) # Time step
end_time = 0.5 # Finish time
time = np.arange(0, end_time, dt) # time range

        ### Initial Conditions:
xo = 0.0
vo = 0.0

        ### Response Histories & Amplitudes:
r_eq = w_eq/wn # Frequency ratio

Y_history = [] # Earh quick input history

Xeq_history = [] # Earhquacke absolute response history
Xeq = (Z(w_eq)[1,1]*(k+complex(0,w_eq*c))*Yo) / ( (Z(w_eq)[0,0]*Z(w_eq)[1,1]) - (Z(w_eq)[0,1])**2 )

Xeqa_history = [] # Absorber absolute response history
Xeqa = - (Z(w_eq)[0,1]*(k+complex(0,w_eq*c))*Yo) / ( (Z(w_eq)[0,0]*Z(w_eq)[1,1]) - (Z(w_eq)[0,1])**2 )

Q_history = [] # Earhquacke relative response history
Q = (Z(w_eq)[1,1]*(m*w_eq**2-ka-complex(0,w_eq*ca))*Yo) / ( (Z(w_eq)[0,0]*Z(w_eq)[1,1]) - (Z(w_eq)[0,1])**2 )

        ### Time Responses:
for t in time:
#--------------------------------------------------
    y = Yo * sin(w_eq*t)                                                      
    Y_history.append(y)                                                   
#--------------------------------------------------
    x_eq = abs(Xeq)*sin(w_eq*t-phase(Xeq)) # Story displacement
    Xeq_history.append(x_eq)
    
    x_eqa = abs(Xeqa)*sin(w_eq*t-phase(Xeqa)) # Absorber displacement
    Xeqa_history.append(x_eqa)
#--------------------------------------------------
Xeq_w_history = np.zeros((len(w_list),1)) # Earhquacke absolute response history
Xeqa_w_history = np.zeros((len(w_list),1)) # Earhquacke absolute response history
for j in range(0,len(w_list)):
        w = w_list[j]
        Xeq = (Z(w)[1,1]*(k+complex(0,w*c))*Yo) / ( (Z(w)[0,0]*Z(w)[1,1]) - (Z(w)[0,1])**2 )
        Xeqa = - (Z(w)[0,1]*(k+complex(0,w*c))*Yo) / ( (Z(w)[0,0]*Z(w)[1,1]) - (Z(w)[0,1])**2 )
        Xeq_w_history[j] = abs(Xeq)
        Xeqa_w_history[j] = abs(Xeqa)

#-------------------------------------------------
#for t in time:
 #   q = Q*sin(w_eq*t-th)
  #  Q_history.append(q)

    
#************************************************************************************************

    ## Results and Responses:

        ### Abs response x_eq(t):
fig1 = plt.figure()
rect = fig1.patch

graphEQ = fig1.add_subplot(1,1,1)
graphEQ.plot(time, Xeq_history, linewidth=2)
graphEQ.plot(time, Xeqa_history, linewidth=2)
graphEQ.plot(time, Y_history, linewidth=1.0)
graphEQ.legend(['x_eq(t) [m]', 'x_eqa(t) [m]', 'y(t) [m]'], fontsize=12, loc='lower right')
graphEQ.set_title('Absolute Response', fontsize=25)
graphEQ.set_xlabel('Time [s]', fontsize=20)
graphEQ.grid(True)

fig2 = plt.figure()
rect = fig2.patch
graphEQ2 = fig2.add_subplot(1,1,1)
graphEQ2.plot(w_list, Xeq_w_history, color='orange', linewidth=6)
graphEQ2.plot(w_list, Xeqa_w_history, color='darkblue', linewidth=2)
graphEQ2.set_title('"Floor & Absorber Amplitude"', fontsize=25)
graphEQ2.set_ylabel('X ', fontsize=20)
graphEQ2.set_xlabel('w_eq', fontsize=20)
graphEQ2.legend(['Xeq(t) [m]', 'Xeqa(t) [m]'], fontsize=13, loc='upper right')
graphEQ2.grid(True)


plt.show()
