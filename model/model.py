
from .parameters import *

from .farm import farm
from .greenhouse import greenhouse
from .ext_clima import ext_clima
from .ComputeFluxes import ComputeFluxes

def model(t,z, climate, daynum):

    gh_state = z[0:11]
    fm_state = z[11:]

    external_radiation, external_climate = ext_clima(t, climate, daynum)

    fluxes = ComputeFluxes(t, gh_state,
                             fm_state, external_radiation)

    dgreenhouse = greenhouse(t,gh_state,
                             fm_state, 
                             external_radiation,
                             external_climate,
                             fluxes)

    dfarm       = farm(t,fm_state,external_radiation,fluxes)

    
    return np.concatenate((dgreenhouse,dfarm))

