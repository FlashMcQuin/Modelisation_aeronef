"""
Author : Hugo Quin

"""
import numpy as np
import sisopy31 as siso
import atm_std as atm
from aeronef import *
import control

aero = Aeronef()
aero.cog(aero.c) #define the value of c
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
state_space = c.ss(A,B,C,D)
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
# x = [gamma, alpha, q, theta, z]
"""
--------------------Step responses in short-period mode : alpha and q ------------------------
"""
sys_a, tf_a, sys_q, tf_q = aero.short_period(A, B)

Ya , Ta = c.matlab.step(sys_a , 10 )
Yq , Tq = c.matlab.step(sys_q , 10 )
plt.plot(Ta ,Ya , 'b' , Tq ,Yq , 'r', lw=2)
step_info1 = c.step_info(Ya,Ta)
Tsa = step_info1['SettlingTime']
Tra = step_info1['RiseTime']
Osa = step_info1['Overshoot']
yya = interp1d(Ta,Ya)
plt.plot(Tsa, yya(Tsa),'bs')
plt.text(Tsa, yya(Tsa)+0.2, round(Tsa,2))
step_info2 = c.step_info(Yq,Tq)
Tsq = step_info2['SettlingTime']
Trq = step_info2['RiseTime']
Osq = step_info2['Overshoot']
yyq = interp1d(Tq,Yq)
plt.plot(Tsq, yyq(Tsq),'rs')
plt.text(Tsq, yyq(Tsq)-0.3, Tsq)
print('alpha Setlling time 5%% = %f s'%Tsa)
print('q Setlling time 5%%     = %f s'%Tsq)
plt.title(r' Step response $\alpha/\delta_m$ and $q/\delta_m$')
plt.legend((r'$\alpha/\delta_m$ ', r'$q/\delta_m$'))
plt.xlabel ( 'Time (s)' )
plt.ylabel ( r'$\alpha$ (rad) & $q$ ( rad/s )' )
plt.show()

"""
------------- Step responses for phugoid mode : v and gamma -------------------------------
"""

ph_v_ss, ph_v_tf, ph_gammma_ss, ph_gammma_tf= aero.phugoid_mode(A,B)
Yv , Tv = c.matlab.step(ph_v_ss, 700)
Yg , Tg = c.matlab.step(ph_gammma_ss, 700 )
plt.plot(Tv ,Yv , 'b' , Tg ,Yg , 'r', lw=2)
step_info1 = c.step_info(Yv,Tv)
Tsv = step_info1['SettlingTime']
Trv = step_info1['RiseTime']
Osv = step_info1['Overshoot']
yyv = interp1d(Tv,Yv)
plt.plot(Tsv, yyv(Tsv),'bs')
plt.text(Tsv, yyv(Tsv)-0.2, Tsv)
step_info2 = c.step_info(Yg,Tg)
Tsg = step_info2['SettlingTime']
Trg = step_info2['RiseTime']
Osg = step_info2['Overshoot']
yyg=interp1d(Tg,Yg)
plt.plot(Tsg, yyg(Tsg),'rs')
plt.text(Tsg, yyg(Tsg)-0.2, Tsg)
print('alpha Setlling time 5%% = %f s'%Tsv)
print('q Setlling time 5%%     = %f s'%Tsg)
plt.title(r' Step response $V/\delta_m$ and $\gamma/\delta_m$')
plt.legend((r'$V/\delta_m$ ', r'$\gamma/\delta_m$'))
plt.xlabel ( 'Time (s)' )
plt.ylabel ( r'$V$ (rad) & $\gamma$ ( rad/s )' )
plt.show()

#Perfect Auto-throttle
A_new = A[1:, 1:] # size 5 without V
B_new = B[1:] # size 5 without V

C_gamma= np.array([[1,0,0,0,0]])
C_alpha=np.array([[0,1,0,0,0]])
C_q= np.array([[0,0,1,0,0]])
C_z= np.array([[0,0,0,0,1]])

aero.open_loop(A_new, B_new, C_q, D, mode = "", title = " q")
plt.show()

"""
-------------------- q feedback loop : --------------------------------
"""

#k=aero.find_k(A_new, B_new, C_q, D)
Kr=-0.08528140111595887

Aq, Bq, sys_q, tf_q= aero.correction_open_loop(A_new, B_new, C_q, D, k=Kr)
print("Aq : ", Aq, "Bq : ", Bq)
Y, t = c.matlab.step(sys_q, 10) # 10 seconds
plt.plot(t, Y)
plt.title("Step response of the closed loop (q feedback loop)")
plt.xlabel("Time (s)")
plt.ylabel("Pitch Rotation Speed (rad/s)")
plt.show()
sys_calpha = c.ss(Aq, Bq, C_alpha, D)
"""
block1 = c.feedback(Kr, sys_q)
block_q = c.series(block1,sys_calpha)
Y, t = c.matlab.step(block_q, 10) # 10 seconds
plt.plot(t, Y, label="alpha with kq")
plt.title("Closed loop alpha")
plt.legend()
"""

""" 
--------Washout filter : ------------
"""
# t = 2/w_sp
to = 2/0.2145
#to = 0.75
washout= c.tf([to, 0], [to,1])
sys_FTBO = tf_a
sys_FTBF = c.series(1/Kr,c.feedback(Kr,tf_q),tf_a) 
sys_FTBF_filter = c.series(1/Kr , c.feedback(Kr,c.series(tf_q, washout)),tf_a) 

T = 10
step_1 , t_1 = c.matlab.step(sys_FTBO,T = T)
step_2 , t_2 = c.matlab.step(sys_FTBF,T = T)
step_3 , t_3 = c.matlab.step(sys_FTBF_filter,T = T)

plt.plot(t_1,step_1,label = "Open Loop")
plt.plot(t_2,step_2,label = "Closed Loop")
plt.plot(t_3,step_3,label = "Closed Loop and Washout Filter")

plt.legend()
plt.title("Step response")
plt.xlabel("Time (in s)")
plt.ylabel(r'$\alpha$')
plt.show()

""" 
--------------------- Gamma feedback loop : ---------------------------
"""

aero.open_loop(Aq, Bq, C_gamma, D, mode = "", title = "of the flight path angle")
plt.show()
#kgamma=aero.find_k(Aq, Bq, C_gamma, D)
kgamma =  17.97283 # second tuning

Agamma, Bgamma, sys_gamma, tf_gamma = aero.correction_open_loop(Aq, Bq, C_gamma, D, k=kgamma)
print("Agamma : ", Agamma, "Bgamma : ", Bgamma)

Y, t = c.matlab.step(tf_gamma, 10) # 10 seconds
plt.plot(t, Y)
plt.title("Step response of the flight path angle hold mode corrected with kgamma")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.show()

""" 
------------------------------ z feedback loop : ------------------------------------
"""
aero.open_loop(Aq, Bq, C_gamma, D, mode = "", title = " z")
plt.show()
#kz = aero.find_k(Agamma, Bgamma, C_z, D) #trouvé
kz= 0.00411
Az, Bz, sys_z, tf_z = aero.correction_open_loop(Agamma, Bgamma, C_z, D, k=kz)
print("Az : ", Az, "Bz : ", Bz)

Y, t = c.matlab.step(sys_z, 10) # 10 seconds
plt.plot(t, Y)
plt.title("Closed loop corrected with kz ")
plt.show()

""" 
----------------------------- Staturation on gamma control loop :  --------------------------------
"""
A_sat = Agamma 
B_sat = Bgamma 
C_sat = np.array([[0, 1, 0, 0, 0]])
D_sat = D
syst_sat =  control.ss(A_sat,B_sat,C_sat,D_sat)
syst_tf_sat = control.tf(syst_sat)

aero.saturation(A_sat,B_sat,C_sat,D_sat)

""" 
------------------------ New C.O.G -------------------------
"""
aero.cog(aero.c_new) #change value of c to c_new
A_new_cog, B_new_cog, C_new_cog, D_new_cog = aero.state_space(X)
A_new_cog = A[1:, 1:] # size 5 without V
B_new_cog = B[1:] # size 5 without V
C_new_cog = np.array([[0, 1, 0, 0, 0]])

sys_new_cog =  control.ss(A_new_cog, B_new_cog, C_new_cog, D_new_cog)
sys_tf_new_cog = control.tf(sys_new_cog)
print("A = ",A_new_cog,"B = ", B_new_cog)
print(sys_tf_new_cog)
control.matlab.damp(sys_new_cog)

plt.figure(10)
q_new , t_new = control.matlab.step(sys_tf_new_cog)
plt.plot(t_new,q_new)
plt.title("Step response of the new system model")
plt.xlabel("Time (s)")
plt.ylabel(r'$\alpha$')
plt.show()

siso.sisotool(sys_tf_new_cog)
kalpha= 0.19729

Aalpha, Balpha, sys_alpha, tf_alpha = aero.correction_open_loop(A_new_cog, B_new_cog, C_new_cog, D, k=kalpha)
y,t=c.matlab.step(sys_alpha, 10)
plt.plot(t,y)
plt.title(r"Step response with control feedback on $\alpha")
plt.xlabel("Time (s)")
plt.ylabel(r'$\alpha$')
plt.show()
