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
        t     = f['t'][:]
        r_0   = f['r_0'][:]
        Sigma = f['Sigma'][:]
        H_gap = f['H_gap'][:]

        r = np.array([ r_0_*x  for r_0_ in r_0 ])

    return t, r, Sigma, H_gap

ia, ib = 0, 2
t, r, Sigma, H_gap = load(ia, ib, 0)

t_ss = t[0] + const.yr*np.array([1e0, 1e1, 1e2, 1e3, 1e4])
j_ss = [ np.where(t >= t_)[0][0]  for t_ in t_ss ]

for j in range(len(j_ss)):
    print("%d: t_ss = %.4e [yr],  H_gap = %.2e [cm] = %g a" % (j,
                                                               (t[j_ss[j]] - t[0])/const.yr,
                                                               H_gap[j_ss[j]], H_gap[j_ss[j]]/(a_ini[ia,ib]*const.AU)))



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
figheight = 15.15 / 2.54          ## convert cm to in
mpl.rc('figure', figsize=[figwidth, figheight])

fig, ax = plt.subplots(nrows=2, ncols=1)


## Log scale
ax_ = ax[0]
ax_.set_title(r"$\alpha = %g$, $\beta = %g$" % (alpha[ia], beta[ib]))
for j in range(len(j_ss)):
    ax_.loglog(r[j_ss[j]]/const.AU, Sigma[j_ss[j]], '-')
ax_.text(1e0, 3e-8, r"$1$")
ax_.text(1.3e0, 2e-7, r"$10$")
ax_.text(1.9e0, 1e-6, r"$100$")
ax_.text(3e0, 3e-6, r"$1000$")
ax_.text(5e0, 7e-6, r"$10000$")
ax_.set_xlabel(r"$r$ [AU]")
ax_.set_xlim(a_ini[ia,ib], 10)
ax_.set_ylabel(r"$\Sigma$ [g$/$cm$^2$]")
ax_.set_ylim(1e-8, 3e-4)

## Linear scale
ax_ = ax[1]
for j in range(len(j_ss)):
    ax_.semilogx(r[j_ss[j]]/const.AU, Sigma[j_ss[j]], '-')
ax_.text(9e-1, 0.4e-6, r"$1$")
ax_.text(1.03e0, 1.5e-6, r"$10$")
ax_.text(1.4e0, 3e-6, r"$100$")
ax_.text(2.5e0, 4e-6, r"$1000$")
ax_.text(5.2e0, 5e-6, r"$10000$")
ax_.set_xlabel(r"$r$ [AU]")
ax_.set_xlim(a_ini[ia,ib], 10)
ax_.set_ylabel(r"$\Sigma$ [g$/$cm$^2$]")
ax_.set_ylim(0, 8e-6)
ax_.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))


plt.tight_layout()
plt.savefig('sigma_evol.eps')
