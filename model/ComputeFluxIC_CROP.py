from .tools.HeatExchange import  convection, radiation
from .parameters import *
from .tools.sat_conc import sat_conc
from .tools.functions import day
from .Growth import Growth

def ComputeFluxes(t,gh_state,fm_state):

    # Values being calculated
    T_c =  gh_state[0]
    T_i =  gh_state[1]
    T_p =  gh_state[2]
    T_f =  gh_state[3]
    T_s1 = gh_state[4]
    T_s2 = gh_state[5]
    T_s3 = gh_state[6]
    T_s4 = gh_state[7]
    C_w =  gh_state[8]
    C_c =  gh_state[9]


    T_m     = fm_state[0]
    T_v     = fm_state[1]
    T_vmean = fm_state[2]
    T_vsum  = fm_state[3]
    C_buf   = fm_state[4]
    C_fruit = fm_state[5]
    C_leaf  = fm_state[6]
    C_stem  = fm_state[7]
    R_fruit = fm_state[8]
    R_leaf  = fm_state[9]
    R_stem  = fm_state[10]




    # External weather and dependent internal parameter values
    p_w = C_w*R*T_i/M_w # Partial pressure of water [Pa]
    rho_i = ((atm - p_w)*M_a + p_w*M_w)/(R*T_i) # Internal density of air [kg/m^3]   
    LAI = SLA*C_leaf # Leaf area index


    # Option for printing progress in hours - uncomment if needed
    #print('Hour', hour)



    # Ventilation


 

    ## Lights
    A_m_wool = 0.75*A_m # Area of mat exposed
    A_m_water = 0.25*A_m # assumed 25% saturated

    ## Convection
    (QV_i_m, QP_i_m, Nu_i_m) = convection(d_m, A_m_wool, T_i, T_m, ias, rho_i, c_i, C_w)
    QP_i_m = A_m_water/A_m_wool * QP_i_m # Factored down

    # Convection internal air -> vegetation
    A_v_exp = LAI*A_v
    (QV_i_v, QP_i_v, Nu_i_v) = convection(d_v, A_v_exp, T_i, T_v, ias, rho_i, c_i, C_w)
    HV = Nu_i_v*lam/d_v
    
    # Convection external air -> cover

    
    # Convection internal air -> tray

    ## Far-IR Radiation

    A_vvf=min(LAI*p_v*A_f,p_v*A_f)
    F_c_v = min((1-F_c_f)*LAI,(1-F_c_f)) # Cover to vegetation
    F_m_c = max((1-LAI),0.0) # Mat to cover
    F_m_v = 1-F_m_c # Mat to vegetation

    # Radiation vegetation to mat
    QR_v_m = radiation(eps_v, eps_m, rho_v, rho_m, F_v_m, F_m_v, A_vvf, T_v, T_m)
    
    # Radiation cover to vegetation
    QR_c_v = radiation(eps_ci, eps_v, rho_ci, rho_v, F_c_v, F_v_c, A_c_roof, T_c, T_v)
    
    # Radiation cover to mat
    F_c_m = max((1-F_c_f)*(1-LAI),0) # Cover to mat
    QR_c_m = radiation(eps_ci, eps_m, rho_ci, rho_m, F_c_m, F_m_c, A_c_roof, T_c, T_m)


    # Radiation vegetation to cover
    QR_v_c = radiation(eps_v, eps_ci, rho_v, rho_ci, F_v_c, F_c_v, A_vvf, T_v, T_c)
    

    # Radiation mat to vegetation
    QR_m_v = radiation(eps_m, eps_v, rho_m, rho_v, F_m_v, F_v_m, A_m, T_m, T_v)
    

    # Radiation tray to mat
    QR_p_m = radiation(eps_p, eps_m, rho_p, rho_m, F_p_m, F_m_p, A_p, T_p, T_m)
    
    # Radiation mat to tray
    
    # Radiation tray to vegetation
    QR_p_v = radiation(eps_p, eps_v, rho_p, rho_v, F_p_v, F_v_p, A_p, T_p, T_v)
    
    QR_v_p = radiation(eps_v, eps_p, rho_v, rho_p, F_v_p, F_p_v, A_vvf, T_v, T_p)


    # Radiation mat to cover
    QR_m_c = radiation(eps_m, eps_ci, rho_m, rho_ci, F_m_c, F_c_m, A_m, T_m, T_c)
        
    # Radiation mat to tray
    QR_m_p = radiation(eps_m, eps_p, rho_m, rho_p, F_m_p, F_p_m, A_m, T_m, T_p)
        
    # Cover to sky


    QD_m_p = (A_m*lam_p/l_m)*(T_m-T_p)

    ##      Solar radiation
    # We first define the solar elevation angle that determines that absorption of solar radiation. Notation: r is direct radiation, f is diffuse radiation, whilst VIS and NIR stand for visible and near infra-red respectively.


    # Radiation from artifical lighting
    # Solar radiation incident on the cover


    # Transmitted solar radiation



    # Solar radiation absorbed by the vegetation
    # Area = A_v i.e. planted area
    # factor QS by A_v/A_f



    ## Photosynthesis model - Vanthoor


    
    MC_buf_i, MC_fruit_i, MC_leaf_i, MC_stem_i, MC_i_buf = Growth(t,gh_state,fm_state,ec_state)
    
    # QT_v_i     = fluxes[0]
    # QR_c_v     = fluxes[1]
    # QV_i_v     = fluxes[2]
    # QR_m_v     = fluxes[3]
    # QR_p_v     = fluxes[4]
    # QR_v_c     = fluxes[5]
    # QR_v_m     = fluxes[6]
    # QR_v_p     = fluxes[7]
    # MC_buf_i   = fluxes[8]
    # MC_fruit_i = fluxes[9]
    # MC_leaf_i  = fluxes[10]
    # MC_stem_i  = fluxes[11]
    # MC_i_buf   = fluxes[12]
    # QV_i_m     = fluxes[13]
    # QP_i_m     = fluxes[14]
    # QR_m_c     = fluxes[15]
    # QR_m_p     = fluxes[16]
    # QD_m_p     = fluxes[17]
    # QS_m_NIR   = fluxes[18]
    # QR_c_m     = fluxes[19]
    # QR_p_m     = fluxes[20]

    # QS_tot_rNIR = fluxes[21]
    # QS_tot_rVIS = fluxes[22]
    # QS_tot_fNIR = fluxes[23]
    # QS_tot_fVIS = fluxes[24]
    # QS_int_rNIR = fluxes[25]
    # QS_int_rVIS = fluxes[26]
    # QS_int_fNIR = fluxes[27]
    # QS_int_fVIS = fluxes[28]
    # QV_e_c      = fluxes[29]
    # QR_c_sk     = fluxes[30]
    # MC_i_e      = fluxes[31]


    return np.array([ QT_v_i, QR_c_v, QV_i_v, QR_m_v, QR_p_v, QR_v_c, QR_v_m,QR_v_p,
                     MC_buf_i, MC_fruit_i, MC_leaf_i, MC_stem_i, MC_i_buf,
                        QV_i_m, QP_i_m, QR_m_c, QR_m_p,QD_m_p,QS_m_NIR,QR_c_m,QR_p_m,
                        QS_tot_rNIR,QS_tot_rVIS,QS_tot_fNIR,QS_tot_fVIS,
                        QS_int_rNIR,QS_int_rVIS,QS_int_fNIR,QS_int_fVIS,
                        QV_e_c,QR_c_sk,QV_i_e,MW_i_e,MC_i_e])