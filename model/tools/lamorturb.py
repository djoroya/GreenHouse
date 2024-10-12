from ..parameters import *

def lamorturb(Gr, Re):

    free = Gr < 1e5
    Nu_G = 0.5 * free * Gr**0.25 + 0.13*(1-free)*Gr**0.33

    forced = Re < 2e4
    Nu_R = 0.6*forced*Re**0.5 + 0.032*(1-forced)*Re**0.8  

    x = Nu_G > Nu_R

    Nu = x*Nu_G + (1-x)*Nu_R

    Sh = x*Nu*Le**0.25 + (1-x)*Nu*Le**0.33

    return(Nu, Sh)