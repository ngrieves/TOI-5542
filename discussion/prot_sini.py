#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 14:42:42 2022

@author: nolangrieves
"""

import numpy as np
import uncertainties as u

rstar_val = 1.058  #Rsun 
rstar_err = 0.036 #Rsun

vrot_star_sini_val = 3.03 #km/s 
vrot_star_sini_err = 0.50 #km/s 

rstar = u.ufloat(rstar_val*695700,rstar_err*695700) #Rsun to km 
vrot_star_sini = u.ufloat(vrot_star_sini_val*60*60*24,vrot_star_sini_err*60*60*24 ) #km/s to km/day

prot_sini = 2*np.pi*rstar/vrot_star_sini

print(r'Prot/sini = %.2f +/- %.2f days'%(prot_sini.nominal_value,prot_sini.std_dev))
