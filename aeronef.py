import numpy as np
import control as c
import atm_std as atm
import sisopy31 as siso
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
class Aeronef():
    def __init__(self, M = 1.21, z=1050):
        self.z0=z*0.3048            #convert feet to meters
        self.m =8400                #kg
        self.Cx0=0.0355
        self.Cz_alpha=2.7
        self.Cz_deltam=0.9
        self.delta_m0=0
        self.alpha0=0.006
        self.f=0.608
        self.f_delta=0.9
        self.k=0.34
        self.Cm_q= -0.48
        self.l_ref = 5.24           #m
        self.l_t= self.l_ref*3/2
        self.g0=9.81
        self.c =0.52
        self.c_new=self.f*1.1
        self.S=34
        self.rg=2.65                #m
        self.tau = 0.75             #s
        # dynamic pressure Pa
        self.gamma0 = 0
        self.alpha0 = 0
        self.theta0 = 0
        self.F=self.f*self.l_t
        self.F_delta=self.f_delta*self.l_t
        self.Iyy=self.m*self.rg**2

    def cog(self, c_val):
        #Centers of gravity
        self.G = c_val*self.l_t
        self.dx = c_val*self.l_t-self.f*self.l_t
        self.dy= c_val*self.l_t - self.f_delta*self.l_t

    def compute_alpha_eq(self, X):
        # X = [V, gamma, aplha, Q, tetha, z]
        alpha0 = X[2]
        Q0 = X[3]
        print(Q0)
        alpha_eq = 0
        F_pxeq = 0

        err = -1
        epislon = 1e-5
        while err > epislon or err == -1:

            C_zeq = (self.m * self.g0 - F_pxeq * np.sin(alpha_eq)) / (Q0 * self.S)
            C_xeq = self.Cx0 + self.k * C_zeq**2

            C_xdeltam = 2 * self.k * C_zeq * self.Cz_deltam

            num = C_xeq * np.sin(alpha_eq) + C_zeq * np.cos(alpha_eq)
            den = C_xdeltam * np.sin(alpha_eq) + self.Cz_deltam * np.cos(alpha_eq)
            delta_meq = self.delta_m0 - num * self.dx / (den * (self.dy - self.dx))

            correction = C_zeq / self.Cz_alpha - self.Cz_deltam * delta_meq / self.Cz_alpha
            alpha_eq_new = alpha0 + correction
            F_pxeq = Q0 * self.S * C_xeq / np.cos(alpha_eq)

            err = abs(alpha_eq_new - alpha_eq)
            alpha_eq = alpha_eq_new

        return C_xeq, C_zeq, alpha_eq, delta_meq, F_pxeq

    def state_space(self, X):
        """Compute state space matrixes of the model.

        Args:
            X (np.array([V, gamma, aplha, Q, tetha, z]))

        Returns:
            tuple(A, B, C, D): State space matrixes.
        """
        gamma_eq = 0
        Ft=0
        Q = X[3]
        self.V_eq = X[0]
        # Calculus simplification variables.
        Cx_eq, Cz_eq, self.alpha_eq, delta_eq, F_eq = self.compute_alpha_eq(X)
        print(self.alpha_eq)
        #print(Cx_eq, Cz_eq, alpha_eq, delta_eq, F_eq)
        QSV = Q*self.S/(self.m*self.V_eq)
        QSI = Q*self.S*self.l_ref/self.Iyy
        Cx_alpha = 2*self.k*Cz_eq*self.Cz_alpha
        Cx_deltam = 2*self.k*Cz_eq*self.Cz_deltam
        Cm_alpha = self.dx*(Cx_alpha*np.sin(self.alpha_eq)+self.Cz_alpha*np.cos(self.alpha_eq))/self.l_ref
        Cm_deltam = self.dy*(Cx_deltam*np.sin(self.alpha_eq)+self.Cz_deltam*np.cos(self.alpha_eq))/self.l_ref
        
        Xv = 2*QSV*Cx_eq
        X_alpha = F_eq*np.sin(self.alpha_eq)/(self.m*self.V_eq) + QSV*Cx_alpha
        X_gamma=self.g0*np.cos(gamma_eq)/self.V_eq
        X_deltam = QSV*Cx_deltam
        Xt=-Ft*np.cos(self.alpha_eq)/(self.m*self.V_eq) # =0

    
        self.m_alpha= QSI*Cm_alpha
        self.m_q=QSI*self.l_ref*self.Cm_q/self.V_eq
        self.m_deltam=QSI*Cm_deltam

        Zv=2*QSV*Cz_eq
        self.Z_alpha=F_eq*np.cos(self.alpha_eq)/(self.m*self.V_eq) + QSV*self.Cz_alpha
        Z_gamma = self.g0*np.sin(gamma_eq)/self.V_eq   # gamma_eq = 0 -> Z_gamma = 0
        self.Z_deltam=QSV*self.Cz_deltam
        Zt= Ft*np.sin(delta_eq)/self.m*self.V_eq # =0

        #print(Xv, X_gamma, X_alpha, Zv, Z_alpha, m_alpha, m_q)
        A = np.array([[-Xv, -X_gamma, -X_alpha, 0, 0, 0],
                    [   Zv, 0, self.Z_alpha, 0, 0, 0],
                    [  -Zv, 0, -self.Z_alpha, 1, 0, 0],
                    [   0, 0, self.m_alpha, self.m_q, 0, 0],
                    [   0, 0, 0, 1, 0, 0],
                    [   0, self.V_eq, 0, 0, 0, 0]])
        
        B = np.array([[0, self.Z_deltam, -self.Z_deltam, self.m_deltam, 0, 0]]).T
        C = np.zeros((1,6))
        C[0,2]=1
        D = np.zeros((1,1))

        return A, B, C, D
    
    def short_period(self, A, B):
        Asp = A[2:4 , 2:4]
        Bsp = B[2:4 , 0:1]
        Csp_alpha = np.matrix ( [[1 , 0]] )
        Csp_q = np.matrix ( [[0 , 1]] )
        Dsp = np.matrix ( [[0]])

        sp_alpha_ss = c.ss(Asp,Bsp,Csp_alpha,Dsp)
        c.matlab.damp(sp_alpha_ss)
        print( "\n\nTransfer function alpha / delta_m =" )
        sp_alpha_tf = c.matlab.ss2tf(sp_alpha_ss)
        print(sp_alpha_tf)
        print ( "Static gain of alpha / delta_m = %f "%(c.dcgain(sp_alpha_ss)))

        sp_q_ss = c.ss(Asp,Bsp,Csp_q,Dsp)
        print( "\n\nTransfer function q / delta_m =" )
        sp_q_tf = c.matlab.ss2tf(sp_q_ss)
        print(sp_q_tf)
        print ( "Static gain of q / delta_m = %f "%(c.dcgain(sp_q_ss)))
        return sp_alpha_ss, sp_alpha_tf, sp_q_ss, sp_q_tf
    
    def phugoid_mode(self, A, B):
        Ap=A[ 0 : 2 , 0 : 2 ]
        Bp=B[ 0 : 2 , 0 : 1 ]
        Cpv = np.matrix ( [[1 , 0]] )
        Cpg = np.matrix ( [[0 , 1]] )
        Dp = np.matrix ( [[0]])
        ph_v_ss =  c.ss(Ap ,Bp ,Cpv ,Dp)
        c.matlab.damp(ph_v_ss)
        print( "\n\nTransfer function V / delta_m =" )
        ph_v_tf = c.tf ( ph_v_ss )
        print( ph_v_tf )
        print ( "Static gain of V / delta_m = %f "%(c.dcgain(ph_v_tf)))

        ph_gammma_ss = c.ss(Ap ,Bp ,Cpg ,Dp)
        print( "\n\nTransfer function gamma / delta_m =" )
        ph_gammma_tf = c.matlab.ss2tf(ph_gammma_ss)
        print(ph_gammma_tf)
        print ( "Static gain of q / delta_m = %f "%(c.dcgain(ph_gammma_ss)))
        return ph_v_ss, ph_v_tf, ph_gammma_ss, ph_gammma_tf
    
    def open_loop(self, A, B, C, D, mode, title):
        ss = c.ss(A, B, C, D)
        ri,a,b,xi,w,st = siso.damp(ss)
        print(f"\n----------------- {mode} Mode :    --------------------------------------\n ")
        for i in st :
            print(i)
        tf = siso.ss2tf(ss)
        print(f"Transfer function : \n ", tf)
        Y, t=c.matlab.step(tf, 10)
        plt.plot(t,Y, label= title)
        plt.title("(Open Loop) Step Response "+title)
        plt.xlabel("Time (s)")
        plt.ylabel("(rad)")

    
    def correction_open_loop(self, A, B, C, D, k):
        # k was found with the use of sisotool on the line above
        print("k = ", k)
        Ak = A-k*np.dot(B,C)
        Bk = k*B
        Dk = k*D
        closed_loop_ss = c.ss(Ak, Bk, C, Dk)
        closed_loop_tf = siso.ss2tf(closed_loop_ss)
        print(f"Transfer function : \n ", closed_loop_tf)
        ri,a,b,xi,w,st = siso.damp(closed_loop_ss)
        for i in st : 
            print(i)
        return Ak, Bk, closed_loop_ss, closed_loop_tf

    def max_incidence(self, A, B, C, D, gamma, alpha_eq,alpha_0):
        delta_nz = 3.1*9.81
        alpha_max = alpha_eq + (alpha_eq - alpha_0)*delta_nz
        syst_ss_alpha = c.ss(A,B*gamma,C,D)
        syst_tf_alpha = c.tf(c.minreal(syst_ss_alpha))
        
        alpha,t_alpha = c.matlab.step(syst_tf_alpha)
        response = max(alpha)-alpha_max
        return response
    
    def saturation(self, A, B, C, D): 
        epsilon = pow(10,-5)
        err = 1
        gamma = 0.1
        i = 0
        while err > epsilon and i<100 :
            gamma_old = gamma
            response_1 = self.max_incidence(A, B, C, D, gamma, self.alpha_eq,self.alpha0)
            response_2 = self.max_incidence(A, B, C, D, gamma+0.01, self.alpha_eq, self.alpha0)
            print(f"response_1 = {response_1} / response_2 = {response_2} ")
            df = (response_2 - response_1)/ 0.01
            print(f"df = {df}")
            gamma = gamma - (response_1/df)
            print(f"gamma = {gamma}")
            err = abs(gamma-gamma_old)
            print(f"Error = {err}")
            i+=1    
        print(f"Gamma max : {(gamma*180/np.pi)}")

