import numpy as np
import sisopy31 as siso
import atm_std as atm
from aeronef import *
import control

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

A,B,C,D = Aeronef().state_space(X)

print(A, "\n")

state_space = control.ss(A,B,C,D)
print(state_space)

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

The system possesses a fast mode at 10.622 rad/s 
with a low damping ratio at 0.215.

The system possesses a slow mode at approx 0.03 rad/s
with a damping ratio at 1.
"""

Aph = A[0:2,0:2]
Bph = B[:2]
C_coupe = C[:2,:2]
D_coupe = D[:2]
Aeronef().transcient_phase(Aph, Bph, C_coupe, D_coupe, mode = "Phugoid")

Asp = A[2:4, 2:4]
Bsp = B[2:4]
Aeronef().transcient_phase(Asp, Bsp, C_coupe, D_coupe, mode = "Short Period")
