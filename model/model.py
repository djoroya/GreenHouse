
from .parameters import *

from .farm import farm
from .greenhouse import greenhouse
from .ext_clima import ext_clima
from .ComputeFluxes import ComputeFluxes

def model(t,z, climate, daynum):

    gh_state = z[0:10]
    fm_state = z[10:]

    ec_state = ext_clima(t, climate, daynum)

    fluxes = ComputeFluxes(t,gh_state,
                             fm_state, 
                             ec_state)

    dgreenhouse = greenhouse(t,gh_state,fluxes)

    dfarm       = farm(t,fm_state,fluxes)

    
    return np.concatenate((dgreenhouse,dfarm))

