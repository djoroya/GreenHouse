
from .parameters import *
from .tools.sat_conc import sat_conc
from .tools.functions import day
from .tools.HeatExchange import convection, radiation
def farm(t,gh_state,fm_state,external_radiation):

    # Values being calculated

    T_c = gh_state[0]
    T_i = gh_state[1]
    T_m = gh_state[2]
    T_p = gh_state[3]
    C_w = gh_state[9]
    C_c = gh_state[10]

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

    # External weather and dependent internal parameter values
    p_w = C_w*R*T_i/M_w # Partial pressure of water [Pa]
    rho_i = ((atm - p_w)*M_a + p_w*M_w)/(R*T_i) # Internal density of air [kg/m^3]   
    LAI = SLA*C_leaf # Leaf area index
    C_c_ppm = C_c*R*T_i/(M_c*atm)*1.e6 # External carbon dioxide concentration [ppm]


    # Option for printing progress of run in days - uncomment out if needed

        
    ## Convection

    
    # Convection internal air -> vegetation
    A_v_exp = LAI*A_v
    (QV_i_v, QP_i_v, Nu_i_v) = convection(d_v, A_v_exp, T_i, T_v, ias, rho_i, c_i, C_w)
    QP_i_v = 0 # No condensation/evaporation - transpiration considered separately
    HV = Nu_i_v*lam/d_v
    

    
    ## Far-IR Radiation

    A_vvf=min(LAI*p_v*A_f,p_v*A_f)
    F_c_v = min((1-F_c_f)*LAI,(1-F_c_f)) # Cover to vegetation
    F_m_c = max((1-LAI),0.0) # Mat to cover
    F_m_v = 1-F_m_c # Mat to vegetation


    
    # Radiation vegetation to cover
    QR_v_c = radiation(eps_v, eps_ci, rho_v, rho_ci, F_v_c, F_c_v, A_vvf, T_v, T_c)
    
    # Radiation vegetation to mat
    QR_v_m = radiation(eps_v, eps_m, rho_v, rho_m, F_v_m, F_m_v, A_vvf, T_v, T_m)
    
    # Radiation vegetation to tray
    QR_v_p = radiation(eps_v, eps_p, rho_v, rho_p, F_v_p, F_p_v, A_vvf, T_v, T_p)
    

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
    QS_al_VIS = 0. 

    # Solar radiation incident on the cover
    QS_tot_rNIR = 0.5*SurfaceArea@external_radiation[0] # Direct 
    QS_tot_rVIS = 0.5*SurfaceArea@external_radiation[0] 
    QS_tot_fNIR = 0.5*SurfaceArea@external_radiation[1]  # Diffuse
    QS_tot_fVIS = 0.5*SurfaceArea@external_radiation[1] 

    # Transmitted solar radiation
    QS_int_rNIR = tau_c_NIR*QS_tot_rNIR # J/s total inside greenhouse
    QS_int_rVIS = tau_c_VIS*QS_tot_rVIS
    QS_int_fNIR = tau_c_NIR*QS_tot_fNIR
    QS_int_fVIS = tau_c_VIS*QS_tot_fVIS 


    # Solar radiation absorbed by the vegetation
    # Area = A_v i.e. planted area
    # factor QS by A_v/A_f

    k_fNIR = 0.27 # Near-IR diffuse extinction coefficient [-]
    a_v_fNIR = 0.65 - 0.65*np.exp(-k_fNIR*LAI) # Near-IR diffuse absorption coefficient [-]

    k_fVIS = 0.85 # Visible diffuse extinction coefficient [-]
    a_v_fVIS = 0.95 - 0.9*np.exp(-k_fVIS*LAI) # Visible diffuse absorption coefficient [-]

    k_rNIR = 0.25 + 0.38*np.exp(-0.12*angle) # Near-IR direct extinction coefficient [-]
    a_v_rNIR = 0.67 - 0.06*np.exp(-0.08*angle) - (0.68 - 0.5*np.exp(-0.11*angle))*np.exp(-k_rNIR*LAI) # Near-IR direct absorption coefficient [-]

    k_rVIS = 0.88 + 2.6*np.exp(-0.18*angle) # Visible direct extinction coefficient [-]
    a_v_rVIS = 0.94 - 0.95*np.exp(-k_rVIS*LAI) # Visible direct absorption coefficient [-]

    QS_v_rNIR = (QS_int_rNIR*(1 - a_obs) + QS_al_NIR)*a_v_rNIR*A_v/A_f
    QS_v_fNIR = (QS_int_fNIR*(1 - a_obs))*a_v_fNIR*A_v/A_f
    QS_v_NIR = (QS_v_rNIR + QS_v_fNIR) # factor as planted area not entire floor

    QS_v_rVIS = (QS_int_rVIS*(1 - a_obs) + QS_al_VIS)*a_v_rVIS*A_v/A_f
    QS_v_fVIS = (QS_int_fVIS*(1 - a_obs))*a_v_fVIS*A_v/A_f
    QS_v_VIS = (QS_v_rVIS + QS_v_fVIS) # Used for photosynthesis calc



    # Solar radiation absorbed by the floor
    # factor by (A_f-A_v)/A_f

    ## Transpiration
    QS_int = (QS_int_rNIR + QS_int_rVIS + QS_int_fNIR + QS_int_fVIS)*(1-a_obs)*A_v/A_f # J/s

    #  Vapour pressure deficit at leaf surface
    xa = C_w/rho_i #[-]
    xv = sat_conc(T_v)/rho_i #[-]
    vpd = atm*(xv/(xv + 0.622) - xa/(xa + 0.622)) # [Pa]

    # Stomatal resistance according to Stanghellini
    x = np.exp(-0.24*LAI) # [-]
    a_v_short = 0.83*(1 - 0.70*x)*(1 + 0.58*x**2)*(0.88 - x**2 + 0.12*x**(8/3)) # [-]Absorption for shortwave radiation
    I_s_bar = QS_int*a_v_short/(2*LAI) # [J/s] Mean radiation interacting with leaf surface

    Heavy_CO2 = I_s_bar > 0.
    r_i_CO2 = 1 + Heavy_CO2*6.1e-7*(C_c_ppm - 200)**2
    Heavy_vpd = vpd/1000 < 0.8
    r_i_vpd = Heavy_vpd*(1 + 4.3*(vpd/1000)**2) + (1 - Heavy_vpd)*3.8
    r_st = 82*((QS_int + 4.3)/(QS_int + 0.54))*(1 + 0.023*(T_v - T_k - 24.5)**2)*r_i_CO2*r_i_vpd #[s/m]

    hL_v_i = 2*LAI*H_fg/(rho_i*c_i)*(Le**(2/3)/HV + r_st/(rho_i*c_i))**(-1)

    QT_St = A_v*hL_v_i*(sat_conc(T_v) - C_w) # J/s

    QT_v_i = max(QT_St,0)
    


    ## Photosynthesis model - Vanthoor

    # Consider photosynthetically active radiation to be visible radiation

    T_25 = T_k + 25. # K

    I_VIS=QS_v_VIS # J/s incident on planted area

    PAR = I_VIS/heat_phot/N_A/A_v 

    # The number of moles of photosynthetically active photons per unit area of planted floor [mol{phot}/m^2/s]
    #J/s/(J/photon)/(photons/mol)/m^2 cf Vanthoor 2.3mumol(photons)/J

    Gamma = max((c_Gamma*(T_v - T_k)/LAI + 20*c_Gamma*(1 - 1/LAI)),0) # The CO2 compensation point [mol{CO2}/mol{air}]
    k_switch = C_buf_max # kg/m^2/s
    h_airbuf_buf = 1/(1 + np.exp(s_airbuf_buf*(C_buf - k_switch)))

    C_c_molar=(C_c/rho_i)*(M_a/M_c)
    C_stom = eta*C_c_molar # Stomatal CO2 concentration [mol{CO2}/mol{air}] 

    J_pot = LAI*J_max_25*np.exp(E_j*(T_v - T_25)/(R*T_v*T_25))*(1 + np.exp((S*T_25 - HH)/(R*T_25)))/(1 + np.exp((S*T_v - HH)/(R*T_v))) # [mol{e}/m^2{floor}s]
    J = (J_pot + alph*PAR - ((J_pot + alph*PAR)**2 - 4*theta*J_pot*alph*PAR)**0.5)/(2*theta)
    P = J*(C_stom - Gamma)/(4*(C_stom + 2*Gamma)) # Photosynthesis rate [mol{CO2}/s]
    Resp = P*Gamma/C_stom # Photorespiration rate

    MC_i_buf = (M_carb*h_airbuf_buf*(P - Resp)) # The net photosynthesis rate [kg{CH2O}/m^2/s]

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
    MC_buf_i = c_fruit_g*MC_buf_fruit + c_leaf_g*MC_buf_leaf + c_stem_g*MC_buf_stem
    #MC_buf_i = 0

    # Maintenance respiration
    MC_fruit_i = (c_fruit_m*Q_10**(0.1*(T_vmean - T_25))*C_fruit*(1 - np.exp(-c_RGR*R_fruit)))
    MC_leaf_i = (c_leaf_m*Q_10**(0.1*(T_vmean - T_25))*C_leaf*(1 - np.exp(-c_RGR*R_leaf)))
    MC_stem_i = (c_stem_m*Q_10**(0.1*(T_vmean - T_25))*C_stem*(1 - np.exp(-c_RGR*R_stem)))

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
