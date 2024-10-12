import numpy as np
from ..parameters import *
def sat_conc(T):

    TC = T - T_k   
    spec_hum = np.exp(11.56 - 4030/(TC + 235))
    air_dens = -0.0046*TC + 1.2978
    a = spec_hum*air_dens
    
    return a