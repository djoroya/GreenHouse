
from .parameters import *
from .tools.sat_conc import sat_conc
from .tools.HeatExchange import convection, radiation, conduction
from .tools.functions import day

def greenhouse(t,gh_state,fm_state,external_radiation,external_climate,fluxes):

    # Values being calculated

    T_c =  gh_state[0]
    T_i =  gh_state[1]
    T_m =  gh_state[2]
    T_p =  gh_state[3]
    T_f =  gh_state[4]
    T_s1 = gh_state[5]
    T_s2 = gh_state[6]
    T_s3 = gh_state[7]
    T_s4 = gh_state[8]
    C_w =  gh_state[9]
    C_c =  gh_state[10]

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

    #
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

    T_ext = external_climate[0]# External air temperature (K)
    T_sk  = external_climate[1]# External sky temperature (K)
    RH_e  = external_climate[3]# External relative humidity
    wind_speed = external_climate[2]# External wind speed (m/s)


    Cw_ext = RH_e * sat_conc(T_ext) # External air moisture content
    p_w = C_w*R*T_i/M_w # Partial pressure of water [Pa]
    rho_i = ((atm - p_w)*M_a + p_w*M_w)/(R*T_i) # Internal density of air [kg/m^3]   
    LAI = SLA*C_leaf # Leaf area index
    C_ce = 4.0e-4*M_c*atm/(R*T_ext) # External carbon dioxide concentration [kg/m^3]


    # Option for printing progress of run in days - uncomment out if needed

    hour = np.floor(t/3600) + 1



    ## Convection
    # Convection external air -> cover

    (QV_e_c, QP_e_c, Nu_e_c ) = convection(d_c, A_c, T_ext, T_c, wind_speed, rho_i, c_i, C_w)

    # Convection internal air -> cover

    (QV_i_c, QP_i_c, Nu_i_c) = convection(d_c, A_c, T_i, T_c, ias, rho_i, c_i, C_w)
    QP_i_c = max(QP_i_c,0) # assumed no evaporation from the cover, only condensation

    # Convection internal air -> floor
    
    (QV_i_f, QP_i_f, Nu_i_f) = convection(d_f, A_f, T_i, T_f, ias, rho_i, c_i, C_w)
    QP_i_f = max(QP_i_f,0) # assumed no evaporation from the floor, only condensation
    
    # Convection internal air -> vegetation
    
    # Convection internal air -> mat
    A_m_wool = 0.75*A_m # Area of mat exposed
    A_m_water = 0.25*A_m # assumed 25% saturated

    (QV_i_m, QP_i_m, Nu_i_m) = convection(d_m, A_m_wool, T_i, T_m, ias, rho_i, c_i, C_w)
    
    QP_i_m = A_m_water/A_m_wool * QP_i_m # Factored down
    
    # Convection internal air -> tray
    (QV_i_p, QP_i_p, Nu_i_p) = convection(d_p, A_p, T_i, T_p, ias, rho_i, c_i, C_w)
    QP_i_p = 0 # Assumed no condensation/evaporation from tray

    ## Far-IR Radiation

    F_c_m = max((1-F_c_f)*(1-LAI),0) # Cover to mat
    F_m_c = max((1-LAI),0.0) # Mat to cover

    # Cover to sky
    QR_c_sk = radiation(eps_ce, 1, 0, 0, 1, 0, A_c, T_c, T_sk)

    # Radiation cover to floor
    QR_c_f = radiation(eps_ci, eps_s, rho_ci, rho_s, F_c_f, F_f_c, A_c_roof, T_c, T_f)

    
    # Radiation cover to mat
    QR_c_m = radiation(eps_ci, eps_m, rho_ci, rho_m, F_c_m, F_m_c, A_c_roof, T_c, T_m)
    
    
    # Radiation mat to cover
    QR_m_c = radiation(eps_m, eps_ci, rho_m, rho_ci, F_m_c, F_c_m, A_m, T_m, T_c)
        
    # Radiation mat to tray
    QR_m_p = radiation(eps_m, eps_p, rho_m, rho_p, F_m_p, F_p_m, A_m, T_m, T_p)
        
    # Radiation tray to mat
    QR_p_m = radiation(eps_p, eps_m, rho_p, rho_m, F_p_m, F_m_p, A_p, T_p, T_m)
    
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
    QD_m_p = (A_m*lam_p/l_m)*(T_m-T_p)

    ## Ventilation
    # Leakage (equations for orifice flow from Awbi, Ventilation of Buildings, Chapter 3)
    wind_speed_H = wind_speed*c*H**a # Wind speed at height H
    wind_pressure = Cp*0.5*rho_i*wind_speed_H**2 # Equals DeltaP for wind pressure
    stack_pressure_diff = rho_i*g*H*(T_i - T_ext)/T_i # DeltaP for stack pressure
 
    Qw = Cd*crack_area*(2*wind_pressure/rho_i)**0.5 # Flow rate due to wind pressure
    Qs = Cd*crack_area*(2*abs(stack_pressure_diff)/rho_i)**0.5 # Flow rate due to stack pressure
    Qt = (Qw**2 + Qs**2)**0.5 # Total flow rate
 
    total_air_flow = Qt*crack_length_total/crack_length 
    R_a_min = total_air_flow/V 

    # Ventilation
    DeltaT_vent = T_i - T_sp_vent
    comp_dtv_low = DeltaT_vent > 0 and DeltaT_vent < 4
    comp_dtv_high = DeltaT_vent >= 4
    R_a = R_a_min + comp_dtv_low*(R_a_max - R_a_min)/4*DeltaT_vent + comp_dtv_high*(R_a_max-R_a_min)

    QV_i_e = R_a*V*rho_i*c_i*(T_i - T_ext) # Internal air to outside air [J/s]
 
    MW_i_e = R_a*(C_w - Cw_ext)

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
    QS_tot_rVIS = 0.5*SurfaceArea@external_radiation[0]
    QS_tot_fNIR = 0.5*SurfaceArea@external_radiation[1] # Diffuse
    QS_tot_fVIS = 0.5*SurfaceArea@external_radiation[1]

    # Transmitted solar radiation
    QS_int_rNIR = tau_c_NIR*QS_tot_rNIR # J/s total inside greenhouse
    QS_int_rVIS = tau_c_VIS*QS_tot_rVIS
    QS_int_fNIR = tau_c_NIR*QS_tot_fNIR
    QS_int_fVIS = tau_c_VIS*QS_tot_fVIS 

    # Solar radiation absorbed by the cover and the obstructions
    QS_c = alph_c*(QS_tot_rNIR + QS_tot_rVIS + QS_tot_fNIR + QS_tot_fVIS) # J/s
    QS_i = a_obs*(QS_int_rNIR + QS_int_rVIS + QS_int_fNIR + QS_int_fVIS) 


    # Solar radiation absorbed by the mat
    a_m_fNIR = 0.05 + 0.91*np.exp(-0.5*LAI) # Near-IR diffuse absorption coefficient [-]
    a_m_rNIR = 0.05 + 0.06*np.exp(-0.08*angle) + (0.92 - 0.53*np.exp(-0.18*angle))*np.exp(-(0.48 + 0.54*np.exp(-0.13*angle))*LAI) # Near-IR direct absorption coefficient [-]

    QS_m_rNIR = (QS_int_rNIR*(1 - a_obs) + QS_al_NIR)*a_m_rNIR*A_v/A_f
    QS_m_fNIR = QS_int_fNIR*(1 - a_obs)*a_m_fNIR*A_v/A_f # W
    QS_m_NIR = (QS_m_rNIR + QS_m_fNIR)


    # Solar radiation absorbed by the floor
    # factor by (A_f-A_v)/A_f

    QS_s_rNIR = QS_int_rNIR*(1-a_obs)*alphS_s*(A_f-A_v)/A_f
    QS_s_fNIR = QS_int_fNIR*(1-a_obs)*alphS_s*(A_f-A_v)/A_f
    QS_s_NIR = QS_s_rNIR + QS_s_fNIR




    ## Dehumidification
    MW_cc_i = 0 # No dehumidification included

    # CO2 exchange with outside
    MC_i_e = (R_a*(C_c - C_ce)) # [kg/m^3/s]

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

    dT_m_dt = (1/(A_m*c_m))*(QV_i_m + QP_i_m - QR_m_v - QR_m_c - QR_m_p - QD_m_p + QS_m_NIR)
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



    dgreenhouse = np.array([dT_c_dt,dT_i_dt,dT_m_dt,dT_p_dt,dT_f_dt,
                     dT_s1_dt,
                     dT_s2_dt,dT_s3_dt,dT_s4_dt,
                     dC_w_dt,
                     dC_c_dt])
    
    return dgreenhouse

