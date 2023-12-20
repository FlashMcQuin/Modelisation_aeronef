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
#fig, axs = plt.subplots(2, 2)

# x = [gamma, alpha, q, theta, z]

#Perfect Auto-throttle
A_new = A[1:, 1:] # size 5 without V
B_new = B[1:] # size 5 without V
print("B_new size : ", B_new.size)
C_gamma= np.array([1,0,0,0,0])
C_alpha=np.array([0,1,0,0,0])
C_q= np.array([0,0,1,0,0])
C_z= np.array([0,0,0,0,1])


aero.open_loop(A_new, B_new, C_alpha, D, mode = "", title = " alpha")


#aero.open_loop(A_new, B_new, C_q, D, mode = "", title= " C_q")



#k=aero.find_k(A_new, B_new, C_q, D)
k=-0.08528140111595887
Aq, Bq = aero.correction_open_loop(A_new, B_new, C_q, D, title = "Calpha", k=k)

sys_cq = c.ss(A_new, B_new, C_q, D) # Aq et Bq
sys_calpha = c.ss(A_new, B_new, C_alpha, D)

block1 = c.feedback(k, sys_cq)
block2 = c.series(block1,sys_calpha)
Y, t = c.matlab.step(block2, 10) # 10 seconds
plt.plot(t, Y, label="alpha with kq")
plt.title("Closed loop alpha")
plt.legend()

############### Washout filter part #######################
to = 0.75
washout= c.tf([to, 0], [to,1])
block1 = c.series(washout,c.ss2tf(sys_cq))
block2 = c.feedback(k, block1)
block3 = c.series(1/k, block2, sys_calpha)
Y, t=c.matlab.step(block3, 10) #10 seconds
plt.plot(t, Y, label="alpha with kq and washout filter")
plt.title("Washout filter")
plt.xlabel("t(s)")
plt.ylabel("alpha (rad)")
plt.legend()
plt.show()
#########################################################


#kgamma=aero.find_k(Aq, Bq, C_gamma, D)
kgamma =  17.97283

Agamma, Bgamma = aero.correction_open_loop(Aq, Bq, C_gamma, D, title = "Cgamma", k=kgamma)
sys_cgamma = c.ss(Aq, Bq, C_gamma, D)

Y, t = c.matlab.step(sys_cgamma, 10) # 10 seconds
plt.plot(t, Y, label="z with kgamma")
plt.title("Closed loop ")
plt.legend()
plt.show()

sys_gamma = c.ss(Aq, Bq, C_gamma, D)
block4 = c.feedback(kgamma, sys_gamma)
Y, t = c.matlab.step(block4, 10) # 10 seconds
plt.plot(t, Y, label="z with kgamma")
plt.title("Closed loop corrected with kgamma ")
plt.legend()
plt.show()


#kz = aero.find_k(Agamma, Bgamma, C_z, D) #trouvé
kz= 0.00392

Az, Bz = aero.correction_open_loop(Agamma, Bgamma, C_z, D, title = "Cgamma", k=kz)

sys_z = c.ss(Az, Bz, C_gamma, D)
block5 = c.feedback(kz, sys_z)
Y, t = c.matlab.step(block5, 10) # 10 seconds
plt.plot(t, Y)
plt.title("Closed loop corrected with kz ")
plt.show()


