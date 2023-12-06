import numpy as np
import control as c
import atm_std as atm
import sisopy31 as siso
import matplotlib.pyplot as plt
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
        self.S=34
        self.rg=2.65                #m
        self.tau = 0.75             #s

        #Centers of gravity
        self.G = self.c*self.l_t
        self.F=self.f*self.l_t
        self.F_delta=self.f_delta*self.l_t

        self.Iyy=self.m*self.rg**2

        self.dx = self.c*self.l_t-self.f*self.l_t
        self.dy= self.c*self.l_t - self.f_delta*self.l_t
        # dynamique pressure Pa
        self.gamma0 = 0
        self.alpha0 = 0
        self.theta0 = 0

    def alpha_eq(self, X):
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
        V_eq = X[0]
        # Calculus simplification variables.
        Cx_eq, Cz_eq, alpha_eq, delta_eq, F_eq = self.alpha_eq(X)
        print(alpha_eq)
        #print(Cx_eq, Cz_eq, alpha_eq, delta_eq, F_eq)
        QSV = Q*self.S/(self.m*V_eq)
        QSI = Q*self.S*self.l_ref/self.Iyy
        Cx_alpha = 2*self.k*Cz_eq*self.Cz_alpha
        Cx_deltam = 2*self.k*Cz_eq*self.Cz_deltam
        Cm_alpha = self.dx*(Cx_alpha*np.sin(alpha_eq)+self.Cz_alpha*np.cos(alpha_eq))/self.l_ref
        Cm_deltam = self.dy*(Cx_deltam*np.sin(alpha_eq)+self.Cz_deltam*np.cos(alpha_eq))/self.l_ref
        
        Xv = 2*QSV*Cx_eq
        X_alpha = F_eq*np.sin(alpha_eq)/(self.m*V_eq) + QSV*Cx_alpha
        X_gamma=self.g0*np.cos(gamma_eq)/V_eq
        X_deltam = QSV*Cx_deltam
        Xt=-Ft*np.cos(alpha_eq)/(self.m*V_eq) # =0

    
        m_alpha= QSI*Cm_alpha
        m_q=QSI*self.l_ref*self.Cm_q/V_eq
        m_deltam=QSI*Cm_deltam

        Zv=2*QSV*Cz_eq
        Z_alpha=F_eq*np.cos(alpha_eq)/(self.m*V_eq) + QSV*self.Cz_alpha
        Z_gamma = self.g0*np.sin(gamma_eq)/V_eq   # gamma_eq = 0 -> Z_gamma = 0
        Z_deltam=QSV*self.Cz_deltam
        Zt= Ft*np.sin(delta_eq)/self.m*V_eq # =0

        #print(Xv, X_gamma, X_alpha, Zv, Z_alpha, m_alpha, m_q)
        A = np.array([[-Xv, -X_gamma, -X_alpha, 0, 0, 0],
                    [   Zv, 0, Z_alpha, 0, 0, 0],
                    [  -Zv, 0, -Z_alpha, 1, 0, 0],
                    [   0, 0, m_alpha, m_q, 0, 0],
                    [   0, 0, 0, 1, 0, 0],
                    [   0, V_eq, 0, 0, 0, 0]])
        
        B = np.array([0, Z_deltam, -Z_deltam, m_deltam, 0, 0])
        B=B.reshape(-1,1)
        C = np.zeros((1,6))
        C[0,2]=1
        D = np.zeros((1,1))

        return A, B, C, D
    
    def transcient_phase_open_loop(self,A, B, C, D, mode, title):
        ss = c.ss(A, B, C, D)
        ri,a,b,xi,w,st = siso.damp(ss)
        print(f"\n----------------- {mode} Mode :    --------------------------------------\n ")
        for i in st : 
            print(i)

        tf = siso.ss2tf(ss)
        print(f"Transfer function : \n ", tf)

        Y, t=c.matlab.step(tf)
        plt.plot(t,Y)
        plt.title("Step Response"+title)
        plt.xlabel("Time (s)")
        plt.ylabel("(rad)")
        plt.show()

    def transient_phase_closed_loop(self, A, B, C, D, title):
        ss = c.ss(A, B, C, D)
        #k = siso.sisotool(ss)
        k = -0.08528140111595887
        print("k = ", k)
        Ak = A-k*B@C
        Bk = k*B
        Dk=k*D
        closed_loop_ss = c.ss(Ak, Bk, C, Dk)
        closed_loop_tf = siso.ss2tf(closed_loop_ss)
        print(f"Transfer function : \n ", closed_loop_tf)
        ri,a,b,xi,w,st = siso.damp(closed_loop_ss)
        for i in st : 
            print(i)

        Y, t=c.matlab.step(closed_loop_tf, 10)
        plt.plot(t,Y)
        plt.title("Step Response "+title)
        plt.xlabel("Time (s)")
        plt.ylabel("y(t)")
        plt.show()


