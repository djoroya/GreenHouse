from .tools.HeatExchange import  convection, radiation
from .parameters import *
from .tools.sat_conc import sat_conc
from .tools.functions import day
from .Growth import Growth

def ComputeFluxes(t,gh_state,fm_state,ec_state):

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


    T_ext   = ec_state[0]# External air temperature (K)
    T_sk    = ec_state[1]# External sky temperature (K)
    wind_sp = ec_state[2]# External wind speed (m/s)
    RH_e    = ec_state[3]# External relative humidity


    # External weather and dependent internal parameter values
    p_w = C_w*R*T_i/M_w # Partial pressure of water [Pa]
    rho_i = ((atm - p_w)*M_a + p_w*M_w)/(R*T_i) # Internal density of air [kg/m^3]   
    LAI = SLA*C_leaf # Leaf area index
    C_c_ppm = C_c*R*T_i/(M_c*atm)*1.e6 # External carbon dioxide concentration [ppm]


    # Option for printing progress in hours - uncomment if needed
    #print('Hour', hour)

    Cw_ext = RH_e * sat_conc(T_ext) # External air moisture content

    wind_sp_H = wind_sp*c*H**a # Wind speed at height H
    wind_pressure = Cp*0.5*rho_i*wind_sp_H**2 # Equals DeltaP for wind pressure
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

    (QV_e_c, QP_e_c, Nu_e_c ) = convection(d_c, A_c, T_ext, T_c, wind_sp, rho_i, c_i, C_w)
    
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
    QR_c_sk = radiation(eps_ce, 1, 0, 0, 1, 0, A_c, T_c, T_sk)


    QD_m_p = (A_m*lam_p/l_m)*(T_m-T_p)

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
    # Solar radiation incident on the cover
    QS_tot_rNIR = 0.5*SurfaceArea@ec_state[4:12] # Direct 
    QS_tot_rVIS = 0.5*SurfaceArea@ec_state[4:12]
    QS_tot_fNIR = 0.5*SurfaceArea@ec_state[12:20] # Diffuse
    QS_tot_fVIS = 0.5*SurfaceArea@ec_state[12:20]

    # Transmitted solar radiation
    QS_int_rNIR = tau_c_NIR*QS_tot_rNIR # J/s total inside greenhouse
    QS_int_rVIS = tau_c_VIS*QS_tot_rVIS
    QS_int_fNIR = tau_c_NIR*QS_tot_fNIR
    QS_int_fVIS = tau_c_VIS*QS_tot_fVIS 


    # Solar radiation absorbed by the vegetation
    # Area = A_v i.e. planted area
    # factor QS by A_v/A_f

    QS_al_NIR = 0. # no artificial lighting

    # Solar radiation absorbed by the mat
    a_m_fNIR = 0.05 + 0.91*np.exp(-0.5*LAI) # Near-IR diffuse absorption coefficient [-]
    a_m_rNIR = 0.05 + 0.06*np.exp(-0.08*angle) + (0.92 - 0.53*np.exp(-0.18*angle))*np.exp(-(0.48 + 0.54*np.exp(-0.13*angle))*LAI) # Near-IR direct absorption coefficient [-]
    QS_m_rNIR = (QS_int_rNIR*(1 - a_obs) + QS_al_NIR)*a_m_rNIR*A_v/A_f
    QS_m_fNIR = QS_int_fNIR*(1 - a_obs)*a_m_fNIR*A_v/A_f # W

    QS_m_NIR = (QS_m_rNIR + QS_m_fNIR)
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

 
    C_ce = 4.0e-4*M_c*atm/(R*T_ext) # External carbon dioxide concentration [kg/m^3]

    MC_i_e = (R_a*(C_c - C_ce)) # [kg/m^3/s]

    
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