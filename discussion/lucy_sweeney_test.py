#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 16:06:37 2022

@author: nolangrieves
"""

import numpy as np

ecc = 0.018
ecc_err = 0.012

lctest = 1 - np.exp(-(ecc**2)/(2*ecc_err**2))
print('Lucy & Sweeney Test statistical significance:',lctest)

# CPL equilibrium rotation period
# Porb = 75.12375
# ecc = 0.018
# Peq = Porb/(1.0 + 9.5*ecc**2.0)
# print(Peq)
# Peq = (2.0/3.0)*Porb
# print(Peq)