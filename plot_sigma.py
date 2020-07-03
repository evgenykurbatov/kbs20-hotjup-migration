# -*- coding: utf-8 -*-

import sys
import numpy as np
from numpy import pi, sqrt, exp, sin, cos, tan, log, log10

import h5py

import matplotlib as mpl
import matplotlib.pyplot as plt

from aux import *
from ini import *



def float_to_str(t):
    e = np.trunc(log10(t/const.yr))
    es = r"10^{%g}" % e

    m = t/const.yr / 10**e
    ms = (r"%g \times" % m) if m > 1 else r""

    return r"$" + ms + es + r"$ [yr]"



##
## Read
##

def load(ia, ib, j):
    with h5py.File('migration_%g_%g_%g_%.4f.h5' % (alpha[ia], beta[ib], t_ini[j], a_ini[ia,ib][j]), 'r') as f:
        x = f['x'][:]
        t[ia,ib]     = f['t'][:]
        a[ia,ib]     = f['a'][:]
        r_0[ia,ib]   = f['r_0'][:]
        j_ss[ia,ib]  = f['j_ss'][:]
        Sigma[ia,ib] = f['Sigma'][:]

        r[ia,ib] = np.array([ r_0_*x  for r_0_ in r_0[ia,ib] ])

t     = {}
a     = {}
r_0   = {}
j_ss  = {}
Sigma = {}
r     = {}

load(0, 0, 0)
load(0, 1, 0)
load(0, 2, 0)
load(1, 0, 0)
load(1, 1, 0)
load(1, 2, 0)



##
## Plot
##

## rc settings (see http://matplotlib.sourceforge.net/users/customizing.html#customizing-matplotlib)
mpl.rc('font', family='serif')
mpl.rc('font', size='8.0')
mpl.rc('text', usetex=True)
mpl.rc('lines', linewidth=0.5)
mpl.rc('axes', linewidth=0.5)
mpl.rc('legend', frameon=False)
mpl.rc('legend', handlelength=2.5)

figwidth = 16.0 / 2.54           ## convert cm to in
figheight = 10.2 / 2.54          ## convert cm to in
mpl.rc('figure', figsize=[figwidth, figheight])

fig, ax = plt.subplots(nrows=2, ncols=3)


ia = 0
for ib in range(len(beta)):
    ax_ = ax[ia,ib]
    j = -1
    ax_.semilogx(r[ia,ib][j]/const.AU, Sigma[ia,ib][j], '-')
    legend = [ mpl.lines.Line2D([], [], color='w', label=(r"$\alpha = %g$" % alpha[ia])), \
               mpl.lines.Line2D([], [], color='w', label=(r"$\beta = %g$" % beta[ib])) ]
    ax_.legend(handles=legend, loc='upper right', handlelength=0)
    ax_.set_xlim(0.04, 3e1)
    ax_.set_ylim(ymax=0.045)

ia = 1
for ib in range(len(beta)):
    ax_ = ax[ia,ib]
    j = -1
    ax_.semilogx(r[ia,ib][j]/const.AU, Sigma[ia,ib][j], '-')
    legend = [ mpl.lines.Line2D([], [], color='w', label=(r"$\alpha = %g$" % alpha[ia])), \
               mpl.lines.Line2D([], [], color='w', label=(r"$\beta = %g$" % beta[ib])) ]
    ax_.legend(handles=legend, loc='upper right', handlelength=0)
    ax_.set_xlim(0.04, 3e1)
    ax_.set_ylim(ymax=0.009)

for ax_ in ax[1,:]:
    ax_.set_xlabel(r"$r$ [AU]")

for ax_ in ax[:,0]:
    ax_.set_ylabel(r"$\Sigma$ [g$/$cm$^2$]")


plt.tight_layout()
plt.savefig('sigma.eps')
