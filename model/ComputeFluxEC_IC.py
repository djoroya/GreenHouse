from .tools.HeatExchange import  convection, radiation
from .parameters import *
from .tools.sat_conc import sat_conc


def ComputeFluxEC_IC(t,gh_state,ec_state):

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
    #
    QS_tot_rNIR = gh_state[10]
    QS_tot_rVIS = gh_state[11]
    QS_tot_fNIR = gh_state[12]
    QS_tot_fVIS = gh_state[13]
    QS_int_rNIR = gh_state[14]
    QS_int_rVIS = gh_state[15]
    QS_int_fNIR = gh_state[16]
    QS_int_fVIS = gh_state[17]



    T_ext   = ec_state[0]# External air temperature (K)
    T_sk    = ec_state[1]# External sky temperature (K)
    wind_sp = ec_state[2]# External wind speed (m/s)
    RH_e    = ec_state[3]# External relative humidity


    # External weather and dependent internal parameter values
    p_w = C_w*R*T_i/M_w # Partial pressure of water [Pa]
    rho_i = ((atm - p_w)*M_a + p_w*M_w)/(R*T_i) # Internal density of air [kg/m^3]   

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


    # Convection internal air -> vegetation

    # Convection external air -> cover

    (QV_e_c, QP_e_c, Nu_e_c ) = convection(d_c, A_c, T_ext, T_c, wind_sp, rho_i, c_i, C_w)
    

    # Cover to sky
    QR_c_sk = radiation(eps_ce, 1, 0, 0, 1, 0, A_c, T_c, T_sk)



    ##      Solar radiation
    # We first define the solar elevation angle that determines that absorption of solar radiation. Notation: r is direct radiation, f is diffuse radiation, whilst VIS and NIR stand for visible and near infra-red respectively.


    # Radiation from artifical lighting
    # Solar radiation incident on the cover


    # Transmitted solar radiation



    # Solar radiation absorbed by the vegetation
    # Area = A_v i.e. planted area
    # factor QS by A_v/A_f

    ## Transpiration

    #  Vapour pressure deficit at leaf surface

    # Stomatal resistance according to Stanghellini




    

    ## Photosynthesis model - Vanthoor

 
    C_ce = 4.0e-4*M_c*atm/(R*T_ext) # External carbon dioxide concentration [kg/m^3]

    MC_i_e = (R_a*(C_c - C_ce)) # [kg/m^3/s]

    
    
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
    # QV_e_c      = fluxes[21]
    # QR_c_sk     = fluxes[22]
    # QV_i_e     = fluxes[23]
    # MW_i_e     = fluxes[24]
    # MC_i_e      = fluxes[25]
    # QS_int_rNIR = fluxes[26]
    # QS_int_fNIR = fluxes[27]


    return QV_e_c,QR_c_sk,QV_i_e,MW_i_e,MC_i_e