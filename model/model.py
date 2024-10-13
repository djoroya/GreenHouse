
from .parameters import *

from .farm import farm
from .greenhouse import greenhouse
from .ext_clima import ext_clima
from .ComputeFluxes import ComputeFluxes

def model(t,z, climate, daynum):

    gh_state = z[0:10]
    fm_state = z[10:]

    ec_state = ext_clima(t, climate, daynum)

    QS_tot_rNIR = 0.5*SurfaceArea@ec_state[4:12] # Direct 
    QS_tot_rVIS = 0.5*SurfaceArea@ec_state[4:12]
    QS_tot_fNIR = 0.5*SurfaceArea@ec_state[12:20] # Diffuse
    QS_tot_fVIS = 0.5*SurfaceArea@ec_state[12:20]
    QS_int_rNIR = tau_c_NIR*QS_tot_rNIR # J/s total inside greenhouse
    QS_int_rVIS = tau_c_VIS*QS_tot_rVIS
    QS_int_fNIR = tau_c_NIR*QS_tot_fNIR
    QS_int_fVIS = tau_c_VIS*QS_tot_fVIS 

    gh_state = np.concatenate((gh_state,
                               [QS_tot_rNIR,   
                                QS_tot_rVIS,
                                QS_tot_fNIR,QS_tot_fVIS,
                                QS_int_rNIR,QS_int_rVIS,
                                QS_int_fNIR,QS_int_fVIS]))



    EC_IC_fluxes, CROP_IC_fluxes = ComputeFluxes(t,
                           gh_state,
                           fm_state, 
                            ec_state)

    dgreenhouse = greenhouse(t,gh_state,EC_IC_fluxes, CROP_IC_fluxes)

    dfarm       = farm(t,fm_state, CROP_IC_fluxes)

    
    return np.concatenate((dgreenhouse,dfarm))

