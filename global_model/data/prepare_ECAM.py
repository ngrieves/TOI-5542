#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 10:54:38 2022

@author: nolangrieves
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy import stats
from astropy.time import Time
from barycorrpy import utc_tdb
import os

#load data
ecam = pd.read_csv('original_data/TIC466206508_LC_cleaned_ap20.dat',sep=",",header=0)

##### CONVERT HJD TO BJD #####
hjd = np.array(ecam['HJD_UTC']) + 2450000
np.savetxt('original_data/intermediate/ECAM_TIC466206508_HJDUTC.txt',hjd,fmt=('%.8f'))
os.system('python convert_times.py original_data/intermediate/ECAM_TIC466206508_HJDUTC.txt hjd mid bjd 20:11:11.62 -- -61:08:07.68 lasilla 161')

#load BJD converted data file
bjd_o = np.array(pd.read_csv('original_data/intermediate/ECAM_TIC466206508_HJDUTC.txt.bjd',skiprows=1,sep=r"\s+",header=None)[0])
flux_o = np.array(ecam['flux'])
fluxerr_o = np.array(ecam['err'])

xshift_o = np.array(ecam['xshift'])
yshift_o = np.array(ecam['yshift'])
am_o = np.array(ecam['AM'])
fwhm_o = np.array(ecam['FWHM'])
sky_o = np.array(ecam['sky'])
exptime_o = np.array(ecam['exptime'])

#remove outliers
bjd = bjd_o[flux_o > 0.99]
flux = flux_o[flux_o > 0.99]
fluxerr = fluxerr_o[flux_o > 0.99]

xshift = xshift_o[flux_o > 0.99]
yshift = yshift_o[flux_o > 0.99] 
am = am_o[flux_o > 0.99] 
fwhm = fwhm_o[flux_o > 0.99] 
sky = sky_o[flux_o > 0.99] 
exptime = exptime_o[flux_o > 0.99] 

#check correlations
from scipy.stats.stats import pearsonr
pvalue_flux_am = pearsonr(flux,am)
print('Flux & AM p-value:',pvalue_flux_am)

pvalue_flux_sky = pearsonr(flux,sky)
print('Flux & Sky p-value:',pvalue_flux_sky)

pvalue_flux_bjd = pearsonr(flux,bjd)
print('Flux & BJD p-value:',pvalue_flux_bjd)

pvalue_flux_xshift = pearsonr(flux,xshift)
print('Flux & xshift p-value:',pvalue_flux_xshift)

pvalue_flux_yshift = pearsonr(flux,yshift)
print('Flux & yshift p-value:',pvalue_flux_yshift)

pvalue_flux_fwhm = pearsonr(flux,fwhm)
print('Flux & fwhm p-value:',pvalue_flux_fwhm)

#plot light curve & check correlations
plt.figure(figsize=(12,6))
plt.errorbar(bjd-bjd[0],flux,yerr=fluxerr,
             fmt='k.',markersize=5,linewidth=1,alpha=0.6,label='ECAM')
plt.xlabel('BJD - '+str(bjd[0]),fontsize=20)
plt.ylabel('Flux',fontsize=20)
plt.title('TOI-5542 ECAM',fontsize=24)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout()
plt.savefig('plots/TOI-5542_ECAM.pdf')

plt.figure(figsize=(12,6))
plt.errorbar(am,flux,yerr=fluxerr,
             fmt='k.',markersize=5,linewidth=1,alpha=0.6,label='ECAM')
plt.xlabel('Airmass',fontsize=20)
plt.ylabel('Flux',fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.title('TOI-5542 ECAM',fontsize=24)
plt.text(1.3,0.995,'Pearson Correlation = %.3f \n p-value= %.3e'%(pvalue_flux_am[0],pvalue_flux_am[1]),fontsize=14)
plt.tight_layout()
#plt.savefig('plots/TOI-5542_ECAM_flux_am.pdf')

plt.figure(figsize=(12,6))
plt.errorbar(sky,flux,yerr=fluxerr,
             fmt='k.',markersize=5,linewidth=1,alpha=0.6,label='ECAM')
plt.xlabel('Sky',fontsize=20)
plt.ylabel('Flux',fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.title('TOI-5542 ECAM',fontsize=24)
plt.text(30,0.995,'Pearson Correlation = %.3f \n p-value= %.3e'%(pvalue_flux_sky[0],pvalue_flux_sky[1]),fontsize=14)
plt.tight_layout()
#plt.savefig('plots/TOI-5542_ECAM_flux_sky.pdf')

plt.figure(figsize=(12,6))
plt.errorbar(fwhm,flux,yerr=fluxerr,
             fmt='k.',markersize=5,linewidth=1,alpha=0.6,label='ECAM')
plt.xlabel('FWHM',fontsize=20)
plt.ylabel('Flux',fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.title('TOI-5542 ECAM',fontsize=24)
plt.text(9,0.995,'Pearson Correlation = %.3f \n p-value= %.3e'%(pvalue_flux_fwhm[0],pvalue_flux_fwhm[1]),fontsize=14)
plt.tight_layout()
#plt.savefig('plots/TOI-5542_ECAM_flux_fwhm.pdf')


plt.figure(figsize=(12,6))
plt.errorbar(xshift,flux,yerr=fluxerr,
             fmt='k.',markersize=5,linewidth=1,alpha=0.6,label='ECAM')
plt.xlabel('xshift',fontsize=20)
plt.ylabel('Flux',fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.title('TOI-5542 ECAM',fontsize=24)
plt.text(-6.5,0.995,'Pearson Correlation = %.3f \n p-value= %.3e'%(pvalue_flux_xshift[0],pvalue_flux_xshift[1]),fontsize=14)
plt.tight_layout()
#plt.savefig('plots/TOI-5542_ECAM_flux_xshift.pdf')


plt.figure(figsize=(12,6))
plt.errorbar(yshift,flux,yerr=fluxerr,
             fmt='k.',markersize=5,linewidth=1,alpha=0.6,label='ECAM')
plt.xlabel('yshift',fontsize=20)
plt.ylabel('Flux',fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.title('TOI-5542 ECAM',fontsize=24)
plt.text(17,0.995,'Pearson Correlation = %.3f \n p-value= %.3e'%(pvalue_flux_yshift[0],pvalue_flux_yshift[1]),fontsize=14)
plt.tight_layout()
#plt.savefig('plots/TOI-5542_ECAM_flux_yshift.pdf')


plt.figure(figsize=(12,6))
plt.errorbar(bjd-bjd[0],flux,yerr=fluxerr,fmt='.',color='black',markersize=5,linewidth=1,alpha=0.6,label='ECAM flux')
plt.xlabel('BJD - '+str(bjd[0]),fontsize=20)
plt.ylabel('variable',fontsize=20)
plt.plot(bjd-bjd[0],sky/np.median(sky),'.',color='blue',alpha=0.6,markersize=5,label='sky')
plt.plot(bjd-bjd[0],am/np.median(am),'.',color='red',alpha=0.6,markersize=5,label='am')
plt.legend()
plt.title('TOI-5542 ECAM',fontsize=24)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

np.savetxt('processed_data/TOI-5542_ECAM.txt',np.column_stack([bjd,flux,fluxerr,xshift,yshift,am,fwhm,sky,exptime]),
            fmt=('%.8f','%.6f','%.6f','%.2f','%.2f','%.2f','%.2f','%.2f','%.2f'))

np.savetxt('processed_data/TOI-5542_ECAM_cds.txt',np.column_stack([bjd,flux,fluxerr,xshift,yshift,am,fwhm,sky,exptime]),
            fmt=('%.8f','%.6f','%.6f','%.2f','%.2f','%.2f','%.2f','%.2f','%.2f'),delimiter='\t')
