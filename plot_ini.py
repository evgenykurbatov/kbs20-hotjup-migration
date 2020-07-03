# -*- coding: utf-8 -*-

import sys
import numpy as np
from numpy import pi, sqrt, exp, sin, cos, tan, log, log10

import matplotlib as mpl
import matplotlib.pyplot as plt

from ini import *



##
## Plot
##

## rc settings (see http://matplotlib.sourceforge.net/users/customizing.html#customizing-matplotlib)
mpl.rc('font', family='serif')
mpl.rc('font', size='8.0')
mpl.rc('text', usetex=True)
mpl.rc('lines', linewidth=1.0)
mpl.rc('axes', linewidth=0.5)
mpl.rc('legend', frameon=False)
mpl.rc('legend', handlelength=3)
mpl.rc('legend', markerscale=0)

figwidth = 8.0 / 2.54           ## convert cm to in
figheight = 8.0 / 2.54          ## convert cm to in
#figheight = 0.82*8.0 / 2.54          ## convert cm to in
mpl.rc('figure', figsize=[figwidth, figheight])

#fig, ax = plt.subplots(nrows=2, ncols=1)
fig, ax = plt.subplots(nrows=1, ncols=1)


ax_ = ax
for ia in range(len(alpha)):
    print(alpha[ia])
    print(r"$\alpha = %g$" % alpha[ia])
    ax_.plot(beta, a_ini[ia,:], ls='', marker='.', label=(r"$\alpha = %g$" % alpha[ia]))
ax_.set_xlabel(r"$\beta$")
ax_.set_ylim(0.3, 0.8)
ax_.set_ylabel(r"$a_\mathrm{ini}$ [AU]")
ax_.legend()


plt.tight_layout()
plt.savefig('ini.eps')
