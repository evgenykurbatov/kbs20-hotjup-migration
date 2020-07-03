# -*- coding: utf-8 -*-

import sys
import subprocess
import numpy as np
from numpy import pi, sqrt, exp, sin, cos, tan, log, log10

from ini import *



for ia in range(len(alpha)):
    for ib in range(len(beta)):
        for i in range(len(t_ini)):
            cmd = "python migration.py %g %g %g %g" % (alpha[ia], beta[ib], t_ini[i], a_ini[ia,ib][i])
            print()
            print(cmd, '-'*20)
            subprocess.Popen(cmd)
