# -*- coding: utf-8 -*-

import numpy as np



alpha = [ 1e-3, 1e-2 ]
beta  = [ 0.0, 1.0, 1.5 ]
t_ini = np.array([ 1e7 ])

a_ini = np.empty((len(alpha), len(beta), len(t_ini)))

## 1e-3, 0
a_ini[0,0] = [ 0.4106 ]
## 1e-3, 1
a_ini[0,1] = [ 0.6690 ]
## 1e-3, 1.5
a_ini[0,2] = [ 0.6739 ]

## 1e-2, 0
a_ini[1,0] = [ 0.5338 ]
## 1e-2, 1
a_ini[1,1] = [ 0.6734 ]
## 1e-2, 1.5
a_ini[1,2] = [ 0.6739 ]
