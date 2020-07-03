# -*- coding: utf-8 -*-

import sys
import numpy as np
from numpy import pi, sqrt, exp, sin, cos, tan, log, log10
import scipy as sp
import scipy.integrate

import h5py

from aux import *



##
## Command line options
##

if len(sys.argv) == 5:
    ## Alpha parameter
    alpha = float(sys.argv[1])
    print("alpha = %g" % alpha)
    ## Turbulent viscosity power index
    beta = float(sys.argv[2])
    print("beta = %g" % beta)
    ## Initial time [yr]
    t_ini = float(sys.argv[3]) * const.yr
    print("t_ini = %g [yr]" % (t_ini/const.yr))
    ## Initial orbit [AU]
    a_ini = float(sys.argv[4]) * const.AU
    print("a_ini = %g [AU]" % (a_ini/const.AU))
else:
    print("Use\n\t%s <alpha> <beta> <t_ini [yr]> <a_ini [AU]>" % sys.argv[0])
    print("Example:\n\t%s 1e-3 1.5 1e9 0.5" % sys.argv[0])
    sys.exit(1)



##
## Parameters connected to the star and planet
##

## Final age
t_fin = t_age

## Final orbit
#a_fin = 0.01*const.AU
a_fin = a_ref
#a_fin = 0.1*const.AU

## Initial planetary mass
M_p_ini = const.M_jup



##
## Parameters connected to the gas
##

## Inner radius of the disk
r_0_ini = ( 1 + (M_p_ini/(3*M_s))**(1/3) ) * a_ini
print("ini: r_0 = %.2e [cm] = %.2e [AU] = %.2e a" % (r_0_ini, r_0_ini/const.AU, r_0_ini/a_ini))

## Keplerian orbital frequency of the planet
Omega_0_ini = Omega_K(a_ini)
print("ini: Omega = %.2e [rad/s] = 2pi/(%.2e [yr])" % (Omega_0_ini, 2*pi/Omega_0_ini/const.yr))
## Mass loss rate at the given orbit
dotM_ini = dotM(t_ini, a_ini)
print("ini: dotM = %.2e [g/s] = %.2e M_sol/[yr]" % (dotM_ini, dotM_ini/(const.M_sol/const.yr)))

##
Sigma_0_ini = dotM_ini/(2*pi*r_0_ini**2*Omega_0_ini)
print("ini: Sigma_0 = %.2e [g/cm^2]" % Sigma_0_ini)

## Disk semi-thickness at its inner radius
H_0_ini = cs(T)/Omega_0_ini
print("ini: H_0 = %.2e [cm] = %.2e [AU] = %.2e a" % (H_0_ini, H_0_ini/const.AU, H_0_ini/a_ini))
## Relative disk semi-thickness
h_ini = H_0_ini/r_0_ini
print("ini: h = %g" % h_ini)



##
## The model
##


## Mass flux
def flux(q, mul, h, xi):
    sigma = q[2:]

    ## Torque
    tau = mul * xi*x_**0.5 / (x_ - xi)**2 * (x_**1.5 - xi**1.5) / (x_**0.5 - xi**0.5)**3

    ## Diffusion
    n = alpha*h**2 * np.asarray(x)**beta

    ## The flux
    f = np.concatenate(([1], \
                        - 3 * ( sqrt(x[1:])*n[1:]*sigma[1:] - sqrt(x[:-1])*n[:-1]*sigma[:-1] ) \
                        / (sqrt(x_[1:-1])*dx_[1:-1])
                        + tau[1:-1] * 0.5*(sigma[1:] + sigma[:-1]) , \
                        [None]))

    ## Right boundary condition
    f[-1] = x_[-2]*f[-2]/x_[-1]

    return f


## R.h.s.
def dotq(t, q):
    M_p   = q[0]
    a     = q[1]
    sigma = q[2:]

    r_0 = ( 1 + (M_p/(3*M_s))**(1/3) ) * a
    xi = a/r_0
    Omega_0 = Omega_K(r_0)

    ## Orbital torque factor (see 1980ApJ...241..425G and 2006RPPh...69..119P)
    C_0 = 2.82

    ##
    ## Auxes

    dotM_p = - dotM(t_ini + t, a)

    dota_int = x / (x - xi)**2 * (x**1.5 - xi**1.5) / (x**0.5 - xi**0.5)**3 * sigma
    dota = sp.integrate.simps(dota_int, x, even='avg')
    dota *= 2*a * C_0/pi * M_p*dotM_p/M_s**2 * sqrt(xi)

    ##
    ## Mass transfer

    mul = 2*C_0/pi * (M_p/M_s)**2
    h = cs(T)/Omega_0 / r_0
    f = flux(q, mul, h, xi)

    dlogdotM_p = dlogdotM(t_ini + t, a)
    Q = dlogdotM_p[0] + (dlogdotM_p[1] - 0.5/a) * dota

    dotsigma = - Omega_0 * (x_[1:]*f[1:] - x_[:-1]*f[:-1]) / (x*dx) - Q*sigma

    return np.concatenate(([dotM_p, dota], dotsigma))


## Event marker
def stop_marker(t, q):
    return q[1] - a_fin
stop_marker.terminal  = True
stop_marker.direction = -1


## Positions of the nodes
#x_ = np.logspace(0, log10(10000*const.AU/r_0_ini), 1601)
x_ = np.logspace(0, log10(1e4*const.AU/r_0_ini), 1001)
print("%g <= x_ <= %g" % (x_[0], x_[-1]))
## Positions of the centers
x = 0.5*(x_[1:] + x_[:-1])
## Grid steps
dx = x_[1:] - x_[:-1]
dx_ = np.concatenate(([None], 0.5*(x_[2:] - x_[:-2]), [None]))

## Jacobian sparsity matrix
jac_sparsity = np.zeros((2+x.size, 2+x.size))
## dM/da
jac_sparsity[0,1] = 1
## da/dsigma
jac_sparsity[1,2:] = 1
## dsigma/dM
jac_sparsity[2:,0] = 1
## dsigma/da
jac_sparsity[2:,1] = 1
## dsigma/dsigma
jac_sparsity[2:,2:] = sp.sparse.spdiags(np.ones((3, x.size)), [-1, 0, 1], x.size, x.size).toarray()

## Computational time grid (it's expandable)
_t = np.empty(1)
## Field's grid (it's expandable)
q = np.empty((1, 2+x.size))

## Times to snapshot
tmp = const.yr * np.array([1e7, 1e8, 1e9, 3e9, 4.1e9])
t_ss = np.concatenate([tmp[tmp > t_ini], [t_fin]])
## Indexes in the 't' array corresponding to the snapshot times
j_ss = np.zeros_like(t_ss, dtype=np.int)

## Initial state
q[0,0]  = M_p_ini
q[0,1]  = a_ini
q[0,2:] = np.zeros_like(x)

## Initial time step
dt_ini = 1/Omega_0_ini
print("ini: dt_ini = %.2e [yr]" % (dt_ini/const.yr))
## Logarithmic time step
dlogt = 0.025

## Solve the model
print("Compute mass transfer...")

## Initial time point
_t[0] = 0
## Index of the current time point
j = 0
## Index of the current snapshot time point
jt_ss = 0
## Flag for saving
to_save = False

while True:
    print("%d: t_ini + %e [yr] = %e [yr]" % (j, _t[j]/const.yr, (t_ini + _t[j])/const.yr))

    ## Get the end time point for the current time step
    _t_next = dt_ini * 10**(dlogt*j)
    ## If it is the time to save?
    if t_ini + _t_next >= t_ss[jt_ss]:
        _t_next = t_ss[jt_ss] - t_ini
        to_save = True

    ## Expand the time grid to the next time point
    _t = np.append(_t, [_t_next])
    ## Expand the field grid to the next time point
    q = np.append(q, [np.empty_like(q[0])], axis=0)

    ## Advance to the next time point
    sol = sp.integrate.solve_ivp(dotq, (_t[j], _t[j+1]), q[j], t_eval=[_t[j], _t[j+1]],
                                 method='BDF', events=stop_marker, dense_output=True,
                                 jac_sparsity=jac_sparsity,
                                 atol=1e-6, rtol=1e-3)

    if sol.status == -1:   ## Error occured
        print("\tERROR: sol.status=%d, '%s'" % (sol.status, sol.message))
        break

    if sol.status == 1:    ## Termination event occured
        ## Set current time to an event time
        _t[j+1] = sol.t_events[0][0]
        q[j+1] = sol.sol(_t[j+1])
        print("\tEvent: t_ini + %e [yr] = %e [yr]" % (_t[j+1]/const.yr, (t_ini + _t[j+1])/const.yr))
        print("\t\ta = %g [AU]" % (q[j+1,1]/const.AU))
        ## Snapshot this state
        j_ss[jt_ss] = j+1
        jt_ss += 1
        break

    q[j+1] = sol.y[:,1]
    print("\ta = %g [AU]" % (q[j+1,1]/const.AU))

    if to_save:
        print("\tSave: t_ss[%d] = %e [yr]" % (jt_ss, t_ss[jt_ss]/const.yr))
        j_ss[jt_ss] = j+1
        jt_ss += 1
        to_save = False

    ## If we finished?
    if t_ini + _t[j+1] >= t_fin:
        print("Finished!")
        break

    ## Prepare to the next time step
    j += 1

t = t_ini + _t
print("... done.")



##
## Finalize
##

M_p   = q[:,0]
a     = q[:,1]
sigma = q[:,2:]

r_0 = ( 1 + (M_p/(3*M_s))**(1/3) ) * a
Omega_0 = Omega_K(r_0)
Sigma_0 = dotM(t, a) / (2*pi*r_0**2*Omega_0)
Sigma = np.multiply(Sigma_0, sigma.T).T

r = np.array([ r_0_*x  for r_0_ in r_0 ])

## Gap width at the semi-height level
Sigma_max = Sigma.max(axis=1)
H_gap = np.zeros_like(t)
for j in range(t.size):
    i = np.where(Sigma[j] >= 0.5*Sigma_max[j])[0][0] - 1
    if i >= 1:
        H_gap[j] = (r[j,i+1] - r[j,i]) / (Sigma[j,i+1] - Sigma[j,i]) * (0.5*Sigma_max[j] - Sigma[j,i]) + r[j,i] - r_0[j]
    #print("\t", j, i, Sigma_max[j], H_gap[j]/a[j])



##
## Plot
##

import matplotlib as mpl
import matplotlib.pyplot as plt

fig, ax = plt.subplots(nrows=2, ncols=2)

dashes = [[2,3], [6,3,2,3], [7,3], []]
color = ['#2a608d', '#e77300', '#4daf4a']

ax_ = ax[0,0]
#for j in range(len(t_ss)):
#    ax_.loglog(x, Sigma[j_ss[j],:], dashes=dashes[j])
ax_.semilogx(x, Sigma[-1,:], dashes=dashes[-1])
ax_.set_xlabel(r"$x$")
ax_.set_xlim(xmax=1e3)
ax_.set_ylabel(r"$\Sigma$ [g$/$cm$^2$]")
#ax_.set_ylim(1e-6, 1e-1)

ax_ = ax[1,0]
ax_.semilogx(r[-1]/const.AU, Sigma[-1,:], dashes=dashes[-1])
#ax_.loglog(r[-1]/const.AU, Sigma[-1,:], dashes=dashes[-1])
ax_.set_xlabel(r"$r$ [AU]")
#ax_.set_xlim(xmax=1e3)
ax_.set_ylabel(r"$\Sigma$ [g$/$cm$^2$]")
#ax_.set_ylim(1e-6, 1e-1)

ax_ = ax[0,1]
ax_.axhline(a_fin/const.AU, c='k', ls=':')
ax_.plot(t/const.yr, a/const.AU)
ax_.set_xlabel(r"$t$ [yr]")
ax_.set_ylabel(r"$a$ [AU]")

ax_ = ax[1,1]
#ax_.semilogx(t/const.yr, r_0/a)
#ax_.set_xlabel(r"$t$ [yr]")
#ax_.set_ylabel(r"$r_0/a$")
ax_.semilogx(t/const.yr, H_gap/a)
ax_.set_xlabel(r"$t$ [yr]")
ax_.set_ylabel(r"$H_\mathrm{gap}/a$")

plt.tight_layout()
plt.show()

#sys.exit(0)




##
## Write
##

with h5py.File('migration_%g_%g_%g_%.4f.h5' % (alpha, beta, t_ini/const.yr, a_ini/const.AU), 'w') as f:
    f.create_dataset('M_s',    data=M_s)    .attrs['comment'] = "Stellar mass [g]"
    f.create_dataset('T',      data=T)      .attrs['comment'] = "Gas temperature [K]"
    f.create_dataset('cs',     data=cs(T))  .attrs['comment'] = "Sound velocity [cm s-1]"
    f.create_dataset('alpha',  data=alpha)  .attrs['comment'] = "Alpha-parameter"
    f.create_dataset('beta',   data=beta)   .attrs['comment'] = "Power index for radial dependence of the viscosity coefficient"
    f.create_dataset('x',      data=x)      .attrs['comment'] = "Radial coordinate grid, in the units of 'r_0'."
    f.create_dataset('t',      data=t)      .attrs['comment'] = "Time grid [s]"
    f.create_dataset('M_p',    data=M_p)    .attrs['comment'] = "Planet mass [g]"
    f.create_dataset('a',      data=a)      .attrs['comment'] = "Planetary orbit radius [cm]"
    f.create_dataset('r_0',    data=r_0)    .attrs['comment'] = "Internal radius of the disk [cm]"
    f.create_dataset('H_gap',  data=H_gap)  .attrs['comment'] = "Gap width at the semi-height level [cm]"
    f.create_dataset('Sigma',  data=Sigma)  .attrs['comment'] = "Surface density at the final time [g cm-2]"
    f.create_dataset('j_ss',   data=j_ss)   .attrs['comment'] = "Indexes in the 't' array corresponding to the snapshot times"
