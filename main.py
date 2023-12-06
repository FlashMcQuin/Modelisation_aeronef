import numpy as np
import sisopy31 as siso
import atm_std as atm
from aeronef import *
import control
aero = Aeronef()
M = 1.21

z0=1050*0.3048   #ft
m =8400          #kg

gamma0 = 0
alpha0 = 0
theta0 = 0
hgeo,rho,a = atm.get_cte_atm(z0)
Veq=M*a
Q = (rho*Veq**2) /2

X = np.array([Veq,
              gamma0,
              alpha0,
              Q,
              theta0,
              z0]).T

A,B,C,D = aero.state_space(X)

state_space = control.ss(A,B,C,D)
print("control minreal : ",control.minreal(state_space))
print("State space : ",state_space)


ri,a,b,xi,w,st = siso.damp(state_space)

for i in st : 
    print(i)

"""damp function output : 
0.000-j0.000  xi=1.000  w=0.000 rad/s
0.000-j0.000  xi=1.000  w=0.000 rad/s

-2.279+j10.375  xi=0.215  w=10.622 rad/s
-2.279-j10.375  xi=0.215  w=10.622 rad/s

-0.022-j0.000  xi=1.000  w=0.022 rad/s
-0.048-j0.000  xi=1.000  w=0.048 rad/s

The system possesses a fast mode (short period) at 10.622 rad/s 
with a low damping ratio at 0.215. -> alpha & q

The system possesses a slow mode at approx 0.03 rad/s
with a damping ratio at 1. -> gamma & V (phugoid)
"""

Aph = A[0:2,0:2]
Bph = B[:2]
Cph = C[0,0:2]
Cph_gamma= np.array([1,0])
Cph_alpha=np.array([0,1])
print("taille A et C : ", Aph.size, Cph_gamma.size)
aero.transcient_phase_open_loop(Aph, Bph, Cph_gamma, D, mode = "Phugoid", title = "C_gamma")
aero.transcient_phase_open_loop(Aph, Bph, Cph_alpha, D, mode = "Phugoid", title = "C_alpha")
Asp = A[2:4, 2:4]
Bsp = B[2:4]
Csp= C[0,2:4]
#aero.transcient_phase_open_loop(Asp, Bsp, Csp, D, mode = "Short Period")

A_sv = A[1:, 1:] 
B_sv = B[1:]
Cgamma = np.zeros((1,5))
Cgamma[0,1]=1
aero.transient_phase_closed_loop(A_sv, B_sv, Cgamma, D, title = "Cgamma")

Cq = np.zeros((1,5))
Cq[0,2] = 1
aero.transient_phase_closed_loop(A_sv, B_sv, Cq, D, title = "Cq")



