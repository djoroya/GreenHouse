from .parameters import *

from .ComputeFluxEC_IC import ComputeFluxEC_IC
from .ComputeFluxIC_CROP import ComputeFluxIC_CROP
def ComputeFluxes(t,gh_state,fm_state,ec_state):

    # Values being calculated

    #

    QS_int_rNIR = gh_state[14]
    QS_int_fNIR = gh_state[16]

    QV_e_c,QR_c_sk,QV_i_e,MW_i_e,MC_i_e = ComputeFluxEC_IC(t,gh_state,ec_state)

    QT_v_i, QR_c_v, QV_i_v, QR_m_v, QR_p_v, QR_v_c, QR_v_m,QR_v_p,MC_buf_i, MC_fruit_i, MC_leaf_i, MC_stem_i, MC_i_buf,QV_i_m, QP_i_m, QR_m_c, QR_m_p,QD_m_p,QS_m_NIR,QR_c_m,QR_p_m = ComputeFluxIC_CROP(t,gh_state,fm_state)


    EC_IC_fluxes = np.array([QV_e_c,
                             QR_c_sk,
                             QV_i_e,
                             MW_i_e,
                             MC_i_e])

    CROP_IC_fluxes = np.array([QT_v_i, 
                               QR_c_v, 
                               QV_i_v,
                                QR_m_v, 
                                QR_p_v, 
                                QR_v_c, 
                                QR_v_m,
                                QR_v_p,
                                MC_buf_i, 
                                MC_fruit_i, 
                                MC_leaf_i, 
                                MC_stem_i, 
                                MC_i_buf,
                                QV_i_m, 
                                QP_i_m, 
                                QR_m_c, 
                                QR_m_p,
                                QD_m_p,
                                QS_m_NIR,
                                QR_c_m,
                                QR_p_m,
                                QS_int_rNIR,
                                QS_int_fNIR])

    return EC_IC_fluxes, CROP_IC_fluxes