#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
General Analysis Equations used in study of TOI-5542
Created on Fri Sep  2 16:06:37 2022
@author: nolangrieves
"""

import numpy as np
import uncertainties as u

###########################
### Lucy & Sweeney Test 
###########################
ecc = 0.018
ecc_err = 0.012

lctest = 1 - np.exp(-(ecc**2)/(2*ecc_err**2))
print('Lucy & Sweeney Test statistical significance:',lctest)

###########################
### Prot/sini 
###########################

rstar_val = 1.058  #Rsun 
rstar_err = 0.036 #Rsun

vrot_star_sini_val = 3.03 #km/s 
vrot_star_sini_err = 0.50 #km/s 

rstar = u.ufloat(rstar_val*695700,rstar_err*695700) #Rsun to km 
vrot_star_sini = u.ufloat(vrot_star_sini_val*60*60*24,vrot_star_sini_err*60*60*24 ) #km/s to km/day

prot_sini = 2*np.pi*rstar/vrot_star_sini

print(r'Prot/sini = %.2f +/- %.2f days'%(prot_sini.nominal_value,prot_sini.std_dev))

###########################
### Rossiter-McLaughlin effect
###########################

vsini = 3.03#*km*second**-1
Rp = 1.009*0.102763#jupiterRad.to('solRad')*solRad
Rstar = 1.058#*solRad
b = 0.419
rm = (Rp/Rstar)**2 * np.sqrt(1-b**2) * vsini
#print(rm.to('m*second**-1'))
print('Rossiter-McLaughlin effect: (m/s)',rm)

###########################
### TSM
###########################
Rp_Earth = 1.008 * 11.2089
Tstar = 5700
Rstar = 1.058 # in solar radii
a = 0.330 * 215.032 # in solar radii
Teq = Tstar * np.sqrt(Rstar/a) * (1/4)**(1/4)
Mp_Earth = 1.33 * 317.907
Jmag = 11.266

ScaleFactor = 1.15


TSM = ScaleFactor * (Rp_Earth**3 * Teq) / (Mp_Earth * Rstar**2) * 10**(-Jmag/5)

print('TSM:',TSM)

###########################
###########################

# CPL equilibrium rotation period
# Porb = 75.12375
# ecc = 0.018
# Peq = Porb/(1.0 + 9.5*ecc**2.0)
# print(Peq)
# Peq = (2.0/3.0)*Porb
# print(Peq)