#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Process raw TESS data
Created on Thu Sep  1 10:49:16 2022
@author: nolangrieves
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits


tess_s13 = fits.open('original_data/hlsp_tess-spoc_tess_phot_0000000466206508-s0013_tess_v1_lc.fits')

time_tess_s13 = tess_s13[1].data['TIME'] + 2457000
time_tess_s13 = time_tess_s13 #- 2450000
flux_tess_s13 = (tess_s13[1].data['PDCSAP_FLUX']/np.nanmedian(tess_s13[1].data['PDCSAP_FLUX']))
fluxerr_tess_s13_original = (tess_s13[1].data['PDCSAP_FLUX_ERR']/np.nanmedian(tess_s13[1].data['PDCSAP_FLUX']))

#remove NANS
time_tess_s13 = time_tess_s13[np.isnan(fluxerr_tess_s13_original) == False]
flux_tess_s13 = flux_tess_s13[np.isnan(fluxerr_tess_s13_original) == False]
fluxerr_tess_s13 = fluxerr_tess_s13_original[np.isnan(fluxerr_tess_s13_original) == False]

plt.figure(figsize=(12,10))
plt.title('TESS Sector 13')
plt.errorbar(time_tess_s13,flux_tess_s13,yerr=fluxerr_tess_s13,fmt='k.',markersize=2,alpha=0.6)
plt.xlabel('BJD',fontsize=20)
plt.ylabel('Flux',fontsize=20)
plt.title('TOI-5542 TESS S13',fontsize=24)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout()
plt.savefig('plots/TOI-5542_TESS_S13.pdf')


tess_s27 = fits.open('original_data/tess2020186164531-s0027-0000000466206508-0189-s_lc.fits')

time_tess_s27 = tess_s27[1].data['TIME'] + 2457000
time_tess_s27 = time_tess_s27 #- 2450000
flux_tess_s27 = (tess_s27[1].data['PDCSAP_FLUX']/np.nanmedian(tess_s27[1].data['PDCSAP_FLUX']))
fluxerr_tess_s27_original = (tess_s27[1].data['PDCSAP_FLUX_ERR']/np.nanmedian(tess_s27[1].data['PDCSAP_FLUX']))

#remove NANS
time_tess_s27 = time_tess_s27[np.isnan(fluxerr_tess_s27_original) == False]
flux_tess_s27 = flux_tess_s27[np.isnan(fluxerr_tess_s27_original) == False]
fluxerr_tess_s27 = fluxerr_tess_s27_original[np.isnan(fluxerr_tess_s27_original) == False]

plt.figure(figsize=(12,10))
plt.title('TESS Sector 27')
plt.errorbar(time_tess_s27,flux_tess_s27,yerr=fluxerr_tess_s27,fmt='k.',markersize=2,alpha=0.6)
plt.xlabel('BJD',fontsize=20)
plt.ylabel('Flux',fontsize=20)
plt.title('TOI-5542 TESS S27',fontsize=24)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout()
plt.savefig('plots/TOI-5542_TESS_S27.pdf')


#SAVE
np.savetxt('processed_data/TOI-5542_TESS_S13.txt',np.column_stack([time_tess_s13,flux_tess_s13,fluxerr_tess_s13]),
            fmt=('%.8f','%.8f','%.8f'))


np.savetxt('processed_data/TOI-5542_TESS_S27.txt',np.column_stack([time_tess_s27,flux_tess_s27,fluxerr_tess_s27]),
            fmt=('%.8f','%.8f','%.8f'))


