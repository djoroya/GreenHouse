from .parameters import *
from .tools.functions import day

def Growth(t,gh_state,fm_state,QS_int_rVIS,QS_int_fVIS):

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
    QS_al_VIS = 0. 

    # Solar radiation incident on the cover



    # Solar radiation absorbed by the vegetation
    # Area = A_v i.e. planted area
    # factor QS by A_v/A_f


    k_fVIS = 0.85 # Visible diffuse extinction coefficient [-]
    a_v_fVIS = 0.95 - 0.9*np.exp(-k_fVIS*LAI) # Visible diffuse absorption coefficient [-]


    k_rVIS = 0.88 + 2.6*np.exp(-0.18*angle) # Visible direct extinction coefficient [-]
    a_v_rVIS = 0.94 - 0.95*np.exp(-k_rVIS*LAI) # Visible direct absorption coefficient [-]


    QS_v_rVIS = (QS_int_rVIS*(1 - a_obs) + QS_al_VIS)*a_v_rVIS*A_v/A_f
    QS_v_fVIS = (QS_int_fVIS*(1 - a_obs))*a_v_fVIS*A_v/A_f
    QS_v_VIS = (QS_v_rVIS + QS_v_fVIS) # Used for photosynthesis calc



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


    # MC_buf_i   = fluxes[1]
    # MC_fruit_i = fluxes[2]
    # MC_leaf_i  = fluxes[3]
    # MC_stem_i  = fluxes[4]
    # MC_i_buf   = fluxes[5]



    return MC_buf_i, MC_fruit_i, MC_leaf_i, MC_stem_i, MC_i_buf