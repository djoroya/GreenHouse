
from .parameters import *
from .tools.HeatExchange import convection, radiation, conduction

def greenhouse(t,gh_state,EC_IC_fluxes, CROP_IC_fluxes):

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

    QS_tot_rNIR = gh_state[10]
    QS_tot_rVIS = gh_state[11]
    QS_tot_fNIR = gh_state[12]
    QS_tot_fVIS = gh_state[13]
    QS_int_rNIR = gh_state[14]
    QS_int_rVIS = gh_state[15]
    QS_int_fNIR = gh_state[16]
    QS_int_fVIS = gh_state[17]

    #
        # Fluxes

    QT_v_i     = CROP_IC_fluxes[0]
    QV_i_v     = CROP_IC_fluxes[2]
    MC_buf_i   = CROP_IC_fluxes[8]
    MC_fruit_i = CROP_IC_fluxes[9]
    MC_leaf_i  = CROP_IC_fluxes[10]
    MC_stem_i  = CROP_IC_fluxes[11]
    MC_i_buf   = CROP_IC_fluxes[12]
    QV_i_m     = CROP_IC_fluxes[13]
    QP_i_m     = CROP_IC_fluxes[14]
    QD_m_p     = CROP_IC_fluxes[17]

    QR_c_v     = CROP_IC_fluxes[1]
    QR_p_v     = CROP_IC_fluxes[4]
    QR_c_m     = CROP_IC_fluxes[19]
    QR_p_m     = CROP_IC_fluxes[20]


    QV_e_c     = EC_IC_fluxes[0]
    QR_c_sk    = EC_IC_fluxes[1]
    QV_i_e     = EC_IC_fluxes[2]
    MW_i_e     = EC_IC_fluxes[3]
    MC_i_e     = EC_IC_fluxes[4]
    
    QR_m_v     = CROP_IC_fluxes[3]
    QR_v_c     = CROP_IC_fluxes[5]
    QR_v_m     = CROP_IC_fluxes[6]
    QR_v_p     = CROP_IC_fluxes[7]
    QR_m_c     = CROP_IC_fluxes[15]
    QR_m_p     = CROP_IC_fluxes[16]
    QS_m_NIR   = CROP_IC_fluxes[18]



    # External weather and dependent internal parameter values


    p_w = C_w*R*T_i/M_w # Partial pressure of water [Pa]
    rho_i = ((atm - p_w)*M_a + p_w*M_w)/(R*T_i) # Internal density of air [kg/m^3]   


    # Option for printing progress of run in days - uncomment out if needed

    hour = np.floor(t/3600) + 1



    ## Convection
    # Convection external air -> cover

    # Convection internal air -> cover

    (QV_i_c, QP_i_c, Nu_i_c) = convection(d_c, A_c, T_i, T_c, ias, rho_i, c_i, C_w)
    QP_i_c = max(QP_i_c,0) # assumed no evaporation from the cover, only condensation

    # Convection internal air -> floor
    
    (QV_i_f, QP_i_f, Nu_i_f) = convection(d_f, A_f, T_i, T_f, ias, rho_i, c_i, C_w)
    QP_i_f = max(QP_i_f,0) # assumed no evaporation from the floor, only condensation
    
    # Convection internal air -> vegetation
    
    # Convection internal air -> mat
    
    
    # Convection internal air -> tray
    (QV_i_p, QP_i_p, Nu_i_p) = convection(d_p, A_p, T_i, T_p, ias, rho_i, c_i, C_w)
    QP_i_p = 0 # Assumed no condensation/evaporation from tray

    ## Far-IR Radiation



    # Radiation cover to floor
    QR_c_f = radiation(eps_ci, eps_s, rho_ci, rho_s, F_c_f, F_f_c, A_c_roof, T_c, T_f)

    
    
    
    
    # Radiation tray to floor
    QR_p_f = radiation(eps_p, eps_s, rho_p, rho_s, F_p_f, F_f_p, A_p, T_p, T_f)
       
    # Radiation floor to cover
    QR_f_c = radiation(eps_s, eps_ci, rho_s, rho_ci, F_f_c, F_c_f, A_f, T_f, T_c)
    
    # Radiation floor to tray
    QR_f_p = radiation(eps_s, eps_p, rho_s, rho_p, F_f_p, F_p_f, A_f, T_f, T_p)
    
    
    ## Conduction
    # Conduction through floor
    QD_sf1 = conduction(A_f, lam_s[0], l_s[0], T_f, T_s1)
    QD_s12 = conduction(A_f, lam_s[1], l_s[1], T_s1, T_s2)
    QD_s23 = conduction(A_f, lam_s[2], l_s[2], T_s2, T_s3)
    QD_s34 = conduction(A_f, lam_s[3], l_s[3], T_s3, T_s4)
    QD_s45 = conduction(A_f, lam_s[4], l_s[4], T_s4, T_ss)
    
    # Conduction mat to tray



    ##      Solar radiation
    # We first define the solar elevation angle that determines that absorption of solar radiation. Notation: r is direct radiation, f is diffuse radiation, whilst VIS and NIR stand for visible and near infra-red respectively.



    # Solar radiation absorbed by the cover and the obstructions
    QS_c = alph_c*(QS_tot_rNIR + QS_tot_rVIS + QS_tot_fNIR + QS_tot_fVIS) # J/s
    QS_i = a_obs*(QS_int_rNIR + QS_int_rVIS + QS_int_fNIR + QS_int_fVIS) 




    # Solar radiation absorbed by the floor
    # factor by (A_f-A_v)/A_f

    QS_s_rNIR = QS_int_rNIR*(1-a_obs)*alphS_s*(A_f-A_v)/A_f
    QS_s_fNIR = QS_int_fNIR*(1-a_obs)*alphS_s*(A_f-A_v)/A_f
    QS_s_NIR = QS_s_rNIR + QS_s_fNIR




    ## Dehumidification
    MW_cc_i = 0 # No dehumidification included

    # CO2 exchange with outside

    day_hour_c=(hour/24-np.floor(hour/24))*24
    track=day_hour_c>6 and day_hour_c<20
    Value=added_CO2/Nz/3600./V

    MC_cc_i=Value*track

    ## Photosynthesis model - Vanthoor

    # Consider photosynthetically active radiation to be visible radiation


    ## ODE equations
    
    # Temperature components
    dT_c_dt = (1/(A_c*cd_c))*(QV_i_c + QP_i_c - QR_c_f - QR_c_v - QR_c_m + QV_e_c - QR_c_sk + QS_c)

    dT_i_dt = (1/(V*rho_i*c_i))*(-QV_i_m - QV_i_v - QV_i_f - QV_i_c - QV_i_e - QV_i_p + QS_i)

    dT_p_dt = (1/(A_p*c_p))*(QD_m_p + QV_i_p + QP_i_p - QR_p_f - QR_p_v - QR_p_m)
    dT_f_dt = (1/(rhod_s[0]*A_f*c_s[0]*l_s[0]))*(QV_i_f + QP_i_f - QR_f_c - QR_f_p - QD_sf1 + QS_s_NIR)
    dT_s1_dt = (1/(rhod_s[1]*c_s[1]*l_s[1]*A_f))*(QD_sf1-QD_s12)
    dT_s2_dt = (1/(rhod_s[2]*c_s[2]*l_s[2]*A_f))*(QD_s12-QD_s23)
    dT_s3_dt = (1/(rhod_s[3]*c_s[3]*l_s[3]*A_f))*(QD_s23-QD_s34)
    dT_s4_dt = (1/(rhod_s[4]*c_s[4]*l_s[4]*A_f))*(QD_s34-QD_s45)

    # Water vapour
    dC_w_dt = (1/(V*H_fg))*(QT_v_i-QP_i_c-QP_i_f-QP_i_m-QP_i_p) - MW_i_e + MW_cc_i
    #dC_wdt = -MW_i_e

    # Carbon Dioxide
    dC_c_dt = MC_cc_i - MC_i_e + (M_c/M_carb)*(A_v/V)*(MC_buf_i + MC_fruit_i + MC_leaf_i + MC_stem_i - MC_i_buf)



    dgreenhouse = np.array([dT_c_dt,dT_i_dt,dT_p_dt,dT_f_dt,
                     dT_s1_dt,
                     dT_s2_dt,dT_s3_dt,dT_s4_dt,
                     dC_w_dt,
                     dC_c_dt])
    
    return dgreenhouse

