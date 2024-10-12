
from .lamorturb import lamorturb
from ..parameters import *
from .sat_conc import sat_conc
def convection(d, A, T1, T2, ias, rho, c, C):
    
    g = 9.81
    nu = 15.1e-6
    lam = 0.025
    
    Gr = (g*d**3)/(T1*nu**2)*abs(T1 - T2)
    Re = ias*d/nu
    (Nu, Sh) = lamorturb(Gr,Re)
    
    QV_1_2 = A*Nu*lam*(T1-T2)/d
    QP_1_2 = A*H_fg/(rho*c)*Sh/Le*lam/d*(C - sat_conc(T2))
    #QP_1_2 = 0
    
    return(QV_1_2, QP_1_2, Nu)

def radiation(eps_1, eps_2, rho_1, rho_2, F_1_2, F_2_1, A_1, T_1, T_2):
    
    sigm = 5.67e-8
    
    k = eps_1*eps_2/(1-rho_1*rho_2*F_1_2*F_2_1)
    QR_1_2 = k*sigm*A_1*F_1_2*(T_1**4 - T_2**4)
    
    return(QR_1_2)

def conduction(A, lam, l, T1, T2):
    QD_12 = (A*lam/l)*(T1-T2)
    
    return(QD_12)