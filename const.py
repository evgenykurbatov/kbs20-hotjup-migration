# -*- coding: utf-8 -*-

#import scipy.constants as spconst


##
## Physical constants
##

## Planck's constant [cm^2/s]
h = 6.626e-27
hbar = 1.055e-27
## Gravitational constant [erg cm/g^2]
G = 6.67e-8
## Speed of light [cm/c]
c = 2.99e10
## Electron charge [(erg cm)^(1/2)]
e = 4.803e-10
## Electron mass [g]
m_e = 9.109e-28
## Proton mass [g]
m_p = 1.673e-24
##Neutron mass [g]
m_n = 1.675e-24
## Hydrogen atom mass
m_H = m_p
## Electron-Volt
EV = 1.6022e-12
eV = EV

## Bohr radius [cm]
a_Bohr = 5.3e-9
## Boltzmann constant [erg/K]
k_B = 1.38e-16
## Avogagro constant [mol^{-1}]
N_Avogadro = 6.02214076e23
## Universal gas constant [erg/mol/K]
## k_B N_Avogadro
R_gas = 8.31446262e7
## Universal gas constant [erg/g/mol/K]
## R_gas/m_H/N_Avogadro
RR_gas = 8.2525343e7

## Stefan-Boltzmann constant [erg/cm^2/s/K^4]
## 2 pi^5 k_B^4 / (15 c^2 h^3)
sigma_SB = 5.67e-5
## Radiation constant (or radiation density constant) [erg/cm^3/K^4]
a_rad = 4.0*sigma_SB/c


##
## Astronomical constants
##

## Solar units
M_sol = 1.99e33
R_sol = 6.96e10
L_sol = 3.827e33
Magn_sol_V   = 4.83
Magn_sol_B   = 5.48
Magn_sol_U   = 5.61
Magn_sol_bol = 4.75

## Planets
M_jup = 9.54e-4 * M_sol
R_jup = 6.99e9
M_earth = 3.00e-6 * M_sol
R_earth = 6.37e8

## Parsec
pc = 3.0857e18
## Astronomical Unit
AU = 1.4960e13
au = AU
## Year
year = 3.1557e7
yr = year

## (erg/s)/(cm^2*Hz)
Jy = 1e-23


##
## Units conversion constants
##

## Length units
nm = 1e-7        ## Nanometers [cm]
Angstrom = 1e-8  ## Angstrom [cm]
## Energy units
Joile = 1e7  ## [erg]
Watt = 1e7   ## [erg/s]
