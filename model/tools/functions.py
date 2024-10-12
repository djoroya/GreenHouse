# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 09:18:07 2021

@author: rmw61
"""
import numpy as np
from ..parameters import *
from .sat_conc import sat_conc

    
def T_ext(t):
    # Weather data

    climate = np.genfromtxt('climate.txt', delimiter=',')
    
    deltaT = 600
    n = int(np.ceil(t/deltaT))
    T_e = climate[n, 0] + T_k
    
    return(T_e)
    

    
def Cw_ext(t):
    # Weather data

    climate = np.genfromtxt('climate.txt', delimiter=',')
    
    deltaT = 600
    n = int(np.ceil(t/deltaT))
    RH_e = climate[n, 1]/100;
    
    Cw_e = RH_e * sat_conc(T_ext(t))
    
    return(Cw_e)

def day(t):
    ## Day
    day_new = np.ceil(t/86400)
    return(day_new)
    