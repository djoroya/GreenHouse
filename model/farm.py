
from .parameters import *
from .tools.functions import day
def farm(t,fm_state,external_radiation,fluxes):

    # Values being calculated



    T_v     = fm_state[0]
    T_vmean = fm_state[1]
    T_vsum  = fm_state[2]
    C_buf   = fm_state[3]
    C_fruit = fm_state[4]
    C_leaf  = fm_state[5]
    C_stem  = fm_state[6]
    R_fruit = fm_state[7]
    R_leaf  = fm_state[8]
    R_stem  = fm_state[9]

            # Fluxes

    QT_v_i     = fluxes[0]
    QR_c_v     = fluxes[1]
    QV_i_v     = fluxes[2]
    QR_m_v     = fluxes[3]
    QR_p_v     = fluxes[4]
    QR_v_c     = fluxes[5]
    QR_v_m     = fluxes[6]
    QR_v_p     = fluxes[7]
    MC_buf_i   = fluxes[8]
    MC_fruit_i = fluxes[9]
    MC_leaf_i  = fluxes[10]
    MC_stem_i  = fluxes[11]
    MC_i_buf   = fluxes[12]

    # External weather and dependent internal parameter values
    LAI = SLA*C_leaf # Leaf area index

    # Option for printing progress of run in days - uncomment out if needed

        
    ## Convection


    ## Far-IR Radiation



    ##      Solar radiation
    # We first define the solar elevation angle that determines that absorption of solar radiation. Notation: r is direct radiation, f is diffuse radiation, whilst VIS and NIR stand for visible and near infra-red respectively.

    gamma = np.deg2rad(360.*(day(t) -  80.)/365.) # Year angle [rad] --- day counts from January 1st
    eqn_time = -7.13*np.cos(gamma) - 1.84*np.sin(gamma) - 0.69*np.cos(2.* gamma) + 9.92*np.sin(2.*gamma) # Equation of time [min]
    az = np.deg2rad(360.*((t/(3600.)%24.) + eqn_time/60. - 12.)/24.) # Azimuth [rad]
    delta = np.deg2rad(0.38 - 0.77*np.cos(gamma) + 23.27*np.cos(gamma)) # Declination angle [rad]
    lat = np.deg2rad(latitude)
    angler = np.arcsin(np.sin(lat)*np.sin(delta) + np.cos(lat)*np.cos(delta)*np.cos(az)) # Angle of elevation [rad]
    angle = np.rad2deg(angler)

    # Radiation from artifical lighting
    QS_al_NIR = 0. # no artificial lighting

    # Solar radiation incident on the cover
    QS_tot_rNIR = 0.5*SurfaceArea@external_radiation[0] # Direct 
    QS_tot_fNIR = 0.5*SurfaceArea@external_radiation[1]  # Diffuse

    # Transmitted solar radiation
    QS_int_rNIR = tau_c_NIR*QS_tot_rNIR # J/s total inside greenhouse
    QS_int_fNIR = tau_c_NIR*QS_tot_fNIR


    # Solar radiation absorbed by the vegetation
    # Area = A_v i.e. planted area
    # factor QS by A_v/A_f

    k_fNIR = 0.27 # Near-IR diffuse extinction coefficient [-]
    a_v_fNIR = 0.65 - 0.65*np.exp(-k_fNIR*LAI) # Near-IR diffuse absorption coefficient [-]


    k_rNIR = 0.25 + 0.38*np.exp(-0.12*angle) # Near-IR direct extinction coefficient [-]
    a_v_rNIR = 0.67 - 0.06*np.exp(-0.08*angle) - (0.68 - 0.5*np.exp(-0.11*angle))*np.exp(-k_rNIR*LAI) # Near-IR direct absorption coefficient [-]


    QS_v_rNIR = (QS_int_rNIR*(1 - a_obs) + QS_al_NIR)*a_v_rNIR*A_v/A_f
    QS_v_fNIR = (QS_int_fNIR*(1 - a_obs))*a_v_fNIR*A_v/A_f
    QS_v_NIR = (QS_v_rNIR + QS_v_fNIR) # factor as planted area not entire floor




    # Solar radiation absorbed by the floor
    # factor by (A_f-A_v)/A_f


    ## Photosynthesis model - Vanthoor

    # Consider photosynthetically active radiation to be visible radiation


    ## Crop growth model

    # Flow of carbohydrates from buffer to fruit, leaves and stem
    C_buf_min = 0.05*C_buf_max
    h_buforg_buf =1/(1 + np.exp(s_buforg_buf*(C_buf - C_buf_min)))

    # inhibition terms need temperatures in oC
    h_T_v = 1/(1 + np.exp(s_min_T*((T_v-T_k) - T_min_v)))/(1 + np.exp(s_max_T*((T_v-T_k) - T_max_v))) 
    h_T_v24 = 1/(1 + np.exp(s_min_T24*((T_vmean-T_k) - T_min_v24)))/(1 + np.exp(s_max_T24*((T_vmean-T_k) - T_max_v24))) 

    h_T_vsum = 0.5*(T_vsum/T_sum_end + ((T_vsum/T_sum_end)**2 + 1e-4)**0.5) - 0.5*(((T_vsum - T_sum_end)/T_sum_end)+(((T_vsum - T_sum_end)/T_sum_end)**2 + 1e-4)**0.5) 

    g_T_v24 = 0.047*(T_vmean - T_k) + 0.06

    MC_buf_fruit = (h_buforg_buf*h_T_v*h_T_v24*h_T_vsum*g_T_v24*rg_fruit) 
    MC_buf_leaf = (h_buforg_buf*h_T_v24*g_T_v24*rg_leaf)
    MC_buf_stem = (h_buforg_buf*h_T_v24*g_T_v24*rg_stem)

    # Growth respiration, which is CO2 leaving the buffer


    C_max_leaf = LAI_max/SLA
    MC_leaf_prune = max(C_leaf - C_max_leaf, 0)

    ## ODE equations
    

    dT_v_dt = (1/(c_v*A_v*msd_v))*(QV_i_v - QR_v_c - QR_v_m - QR_v_p + QS_v_NIR - QT_v_i)


    # Plant growth control
    dT_vmean_dt = 1/86400*(T_v - T_vmean) 
    dT_vsum_dt = 1/86400*(T_v - T_k) 

    # Plant carbon exchange
    dC_buf_dt = MC_i_buf - MC_buf_fruit - MC_buf_leaf - MC_buf_stem - MC_buf_i
    dC_fruit_dt = MC_buf_fruit - MC_fruit_i
    dC_leaf_dt = MC_buf_leaf - MC_leaf_i - MC_leaf_prune
    dC_stem_dt = MC_buf_stem - MC_stem_i

    # Plant growth
    dR_fruit_dt = (dC_fruit_dt/C_fruit - R_fruit)
    dR_leaf_dt = ((dC_leaf_dt + MC_leaf_prune)/C_leaf - R_leaf)
    dR_stem_dt = (dC_stem_dt/C_stem - R_stem)
    
    return np.array([dT_v_dt,dT_vmean_dt,dT_vsum_dt,
                     dC_buf_dt,
                     dC_fruit_dt,dC_leaf_dt,dC_stem_dt,
                     dR_fruit_dt,dR_leaf_dt,dR_stem_dt])
