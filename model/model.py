
from .parameters import *

from .farm import farm
from .greenhouse import greenhouse
from .ext_clima import ext_clima
def model(t,z, climate, daynum):

    gh_state = z[0:11]
    fm_state = z[11:]

    external_radiation, external_climate = ext_clima(t, z, climate, daynum)

    dgreenhouse = greenhouse(t,gh_state,
                             fm_state, 
                             external_radiation,
                             external_climate)

    dfarm       = farm(t,gh_state,
                       fm_state,
                       external_radiation)

    
    return np.concatenate((dgreenhouse,dfarm))

