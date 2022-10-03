#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 16:31:20 2022

@author: nolangrieves
"""

import numpy as np


#TOI-5542b
vsini = 3.03#*km*second**-1
Rp = 1.009*0.102763#jupiterRad.to('solRad')*solRad
Rstar = 1.058#*solRad
b = 0.419
rm = (Rp/Rstar)**2 * np.sqrt(1-b**2) * vsini
#print(rm.to('m*second**-1'))
print('Rossiter-McLaughlin effect: (m/s)',rm)


#TSM


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