import numpy as np
from .parameters import *
from .tools.functions import day
def ext_clima(t,z, climate, daynum):
    

    n = int(np.ceil(t/deltaT)) # count
    daynum.append(day(t)) # Day number

    DirectSolarRad  = climate[n, 4:12]
    DiffusSolarRad = climate[n, 12:20]

    external_radiation = np.array([DirectSolarRad, 
                                   DiffusSolarRad])

    T_ext = climate[n, 0] + T_k # External air temperature (K)
    T_sk  = climate[n, 1] + T_k # External sky temperature (K)
    RH_e   = climate[n, 3]/100 # External relative humidity
    wind_speed = climate[n, 2] # External wind speed (m/s)

    external_climate = np.array([T_ext, T_sk, 
                                 wind_speed, RH_e
                                 ])
    
    return external_radiation, external_climate