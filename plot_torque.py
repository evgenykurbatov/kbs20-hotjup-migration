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
        a[ia,ib]     = f['a'][:]
        Sigma[ia,ib] = f['Sigma'][:]
        H_gap[ia,ib] = f['H_gap'][:]

a     = {}
Sigma = {}
H_gap = {}

#load(0, 0, 0)
load(0, 1, 0)
load(0, 2, 0)
#load(1, 0, 0)
load(1, 1, 0)
load(1, 2, 0)

j = -1
h_gap = [ [ H_gap[ia,ib][j]/a[ia,ib][j]
            for ib in [1,2] ]
          for ia in range(len(alpha)) ]
Sigma_max = [ [ Sigma[ia,ib][j].max()
                for ib in [1,2] ]
              for ia in range(len(alpha)) ]



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

figwidth = 8.0 / 2.54           ## convert cm to in
figheight = 8.2 / 2.54          ## convert cm to in
mpl.rc('figure', figsize=[figwidth, figheight])

fig, ax = plt.subplots(nrows=1, ncols=1)


ax_ = ax
j = -1
for ia, alpha_ in enumerate(alpha):
    ax_.semilogy(h_gap[ia], Sigma_max[ia], ls='', marker='.', label=(r"$\alpha = %g$" % alpha_))
h_iso = np.logspace(log10(0.22), log10(0.64), 50)
Sigma_iso = lambda C : C * (h_iso/0.22)**3
for C in [0.002, 0.006]:
    ax_.semilogy(h_iso, Sigma_iso(C), '-', c='#4daf4a')
ax_.semilogy(h_iso, 0.0063 * (h_iso/0.22)**1.9, '-', c='#f781bf')
ax_.legend()
ax_.set_xlabel(r"$H_\mathrm{gap}/a$")
ax_.set_ylabel(r"$\Sigma$ [g$/$cm$^2$]")


plt.tight_layout()
plt.show()
#plt.savefig('torque.eps')
