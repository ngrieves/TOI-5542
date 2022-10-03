#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 11:05:00 2022

@author: nolangrieves
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy import stats
from astropy.time import Time
from barycorrpy import utc_tdb

msize = 1
lsize = 0.1

direc = 'original_data/'
file = 'TIC-466_20210803_SAAO_BJD.dat'

# The columns are:
# 1. BJD / days (TDB_MID)
# 2. Differential flux
# 3. Error in differential flux
#    (formal, not scaled)
# 4. Target x position / pix
# 5. Target y position / pix
# 6. Airmass
# 7. Target FWHM / arcsecs
# 8. Sky background
# 9. Exposure time / s
#10. JD / days (UTC_MID - 2450000) 

saao = pd.read_csv(direc+file,skiprows=50,sep=r"\s+",header=None)
bjd = np.array(saao[0])
flux = np.array(saao[1])
fluxerr = np.array(saao[2])
xshift = np.array(saao[3])
yshift = np.array(saao[4])
am = np.array(saao[5])
fwhm = np.array(saao[6])
sky = np.array(saao[7])
exptime = np.array(saao[8])

#plot original data
plt.figure(figsize=(12,6))
plt.errorbar(bjd,flux,yerr=fluxerr,
             fmt='k.',markersize=msize,linewidth=lsize,label='SAAO')
plt.xlabel('BJD',fontsize=20)
plt.ylabel('Flux',fontsize=20)
plt.title('TOI-5542 SAAP',fontsize=24)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout()
plt.savefig('plots/TOI-5542_SAAO.pdf')


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


plt.figure(figsize=(12,6))
plt.errorbar(bjd-bjd[0],flux,yerr=fluxerr,
             fmt='k.',markersize=5,linewidth=1,alpha=0.6,label='SAAO')
plt.xlabel('BJD - '+str(bjd[0]),fontsize=20)
plt.ylabel('Flux',fontsize=20)
plt.title('TIC466206508 SAAO',fontsize=24)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout()
#plt.savefig('plots/TIC466206508_SAAO.pdf')

plt.figure(figsize=(12,6))
plt.errorbar(am,flux,yerr=fluxerr,
             fmt='k.',markersize=5,linewidth=1,alpha=0.6,label='SAAO')
plt.xlabel('Airmass',fontsize=20)
plt.ylabel('Flux',fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.title('TIC466206508 SAAO',fontsize=24)
plt.text(1.4,1.01,'Pearson Correlation = %.3f \n p-value= %.3e'%(pvalue_flux_am[0],pvalue_flux_am[1]),fontsize=14)
plt.tight_layout()
#plt.savefig('plots/TIC466206508_SAAO_flux_am.pdf')

plt.figure(figsize=(12,6))
plt.errorbar(sky,flux,yerr=fluxerr,
             fmt='k.',markersize=5,linewidth=1,alpha=0.6,label='SAAO')
plt.xlabel('Sky',fontsize=20)
plt.ylabel('Flux',fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.title('TIC466206508 SAAO',fontsize=24)
plt.text(42,1.01,'Pearson Correlation = %.3f \n p-value= %.3e'%(pvalue_flux_sky[0],pvalue_flux_sky[1]),fontsize=14)
plt.tight_layout()
#plt.savefig('plots/TIC466206508_SAAO_flux_sky.pdf')


plt.figure(figsize=(12,6))
plt.errorbar(fwhm,flux,yerr=fluxerr,
             fmt='k.',markersize=5,linewidth=1,alpha=0.6,label='SAAO')
plt.xlabel('FWHM',fontsize=20)
plt.ylabel('Flux',fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.title('TIC466206508 SAAO',fontsize=24)
plt.text(2.2,1.01,'Pearson Correlation = %.3f \n p-value= %.3e'%(pvalue_flux_fwhm[0],pvalue_flux_fwhm[1]),fontsize=14)
plt.tight_layout()
#plt.savefig('plots/TIC466206508_SAAO_flux_fwhm.pdf')


plt.figure(figsize=(12,6))
plt.errorbar(xshift,flux,yerr=fluxerr,
             fmt='k.',markersize=5,linewidth=1,alpha=0.6,label='SAAO')
plt.xlabel('x position',fontsize=20)
plt.ylabel('Flux',fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.title('TIC466206508 SAAO',fontsize=24)
plt.text(110,0.995,'Pearson Correlation = %.3f \n p-value= %.3e'%(pvalue_flux_xshift[0],pvalue_flux_xshift[1]),fontsize=14)
plt.tight_layout()
#plt.savefig('plots/TIC466206508_SAAO_flux_xshift.pdf')


plt.figure(figsize=(12,6))
plt.errorbar(yshift,flux,yerr=fluxerr,
             fmt='k.',markersize=5,linewidth=1,alpha=0.6,label='SAAO')
plt.xlabel('y position',fontsize=20)
plt.ylabel('Flux',fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.title('TIC466206508 SAAO',fontsize=24)
plt.text(32,0.995,'Pearson Correlation = %.3f \n p-value= %.3e'%(pvalue_flux_yshift[0],pvalue_flux_yshift[1]),fontsize=14)
plt.tight_layout()
#plt.savefig('plots/TIC466206508_SAAO_flux_yshift.pdf')

plt.figure(figsize=(12,8))
plt.errorbar(bjd-bjd[0],flux,yerr=fluxerr,
             fmt='k.',markersize=5,linewidth=1,alpha=0.6,label='SAAO flux')
plt.plot(bjd-bjd[0],0.94+am/40,'.',color='green',alpha=0.6,label=r'airmass (arbitrary scaled and shifted)')
plt.plot(bjd-bjd[0],0.94+am**2/40,'.',color='blue',alpha=0.6,label=r'airmass$^{2}$ (arbitrary scaled and shifted)')
plt.plot(bjd-bjd[0],0.94+sky/1000,'.',color='red',alpha=0.6,label=r'sky (arbitrary scaled and shifted)')
plt.plot(bjd-bjd[0],0.94+fwhm/40,'.',color='orange',alpha=0.6,label=r'FWHM (arbitrary scaled and shifted)')
plt.plot(bjd-bjd[0],0.84+xshift/1000,'.',color='purple',alpha=0.6,label=r'x position (arbitrary scaled and shifted)')
plt.plot(bjd-bjd[0],0.94+yshift/1000,'.',color='cyan',alpha=0.6,label=r'y position (arbitrary scaled and shifted)')
plt.xlabel('BJD - '+str(bjd[0]),fontsize=20)
plt.ylabel('Scaled Variable',fontsize=20)
plt.title('TIC466206508 SAAO',fontsize=24)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=14)
plt.tight_layout()
#plt.savefig('plots/TIC466206508_SAAO_correlations.pdf')

###save data
np.savetxt('processed_data/TOI-5542_SAAO.txt',np.column_stack([(bjd),flux,fluxerr,xshift,yshift,am,fwhm,sky,exptime]),
            fmt=('%.8f','%.8f','%.8f','%.8f','%.8f','%.8f','%.8f','%.8f','%.8f'))

np.savetxt('processed_data/TOI-5542_SAAO_cds.txt',np.column_stack([(bjd),flux,fluxerr,xshift,yshift,am,fwhm,sky,exptime]),
            fmt=('%.8f','%.8f','%.8f','%.8f','%.8f','%.8f','%.8f','%.8f','%.8f'),delimiter='\t')

