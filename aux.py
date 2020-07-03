# -*- coding: utf-8 -*-

import numpy as np
from numpy import pi, sqrt, exp, sin, cos, tan, log, log10

import const



## Mass of the star
M_s = const.M_sol

## Mean molecular weigth (assumed to be constant)
mu = 1

## Atmosphere loss model of Garcia Munoz (2007) [https://ui.adsabs.harvard.edu/abs/2007P&SS...55.1426G]
dotM_ref = 1e12   ## g/s
a_ref    = 0.045 * const.AU
t_age    = 4.6e9 * const.yr
## This power law is applicable for `a > 0.015 AU`:
dotM     = lambda t, a : dotM_ref * t_age/t * (a_ref/a)**2
dlogdotM = lambda t, a : [-1/t, -2/a]

## Keplerian orbital frequency
Omega_K = lambda r : sqrt(const.G*M_s/r**3)

## Temperature of gas
T = 1e4  ## K
## Sound velocity
cs = lambda T : sqrt((5/3)*const.RR_gas/mu*T)
