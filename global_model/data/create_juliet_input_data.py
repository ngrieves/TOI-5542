#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 11:20:32 2022

@author: nolangrieves
"""

import juliet
import matplotlib.pyplot as plt
import pickle
import numpy as np
import pandas as pd


#load TESS Sector 13 data
tess13_data = np.loadtxt('processed_data/TOI-5542_TESS_S13.txt',usecols=[0,1,2],delimiter=' ')

#load TESS Sector 27 data
tess27_data = np.loadtxt('processed_data/TOI-5542_TESS_S27.txt',usecols=[0,1,2],delimiter=' ')

#load ngts data
#ngts_data = np.loadtxt('processed_data/TIC466206508_NGTS.txt',usecols=[0,1,2],delimiter=' ')

#load saao data
saao_data = np.loadtxt('processed_data/TOI-5542_SAAO.txt',usecols=[0,1,2,3,4,5,6,7,8],delimiter=' ')

#load ecam data
ecam_data = np.loadtxt('processed_data/TOI-5542_ECAM.txt',usecols=[0,1,2,3,4,5,6,7,8],delimiter=' ')


#### create dictonaries for input photometry data 
times = {}
fluxes = {}
fluxes_error = {}

times['TESS13'] = tess13_data[:,0] - 2400000.0
times['TESS27'] = tess27_data[:,0] - 2400000.0
#times['NGTS'] = ngts_data[:,0] - 2400000.0
times['SAAO'] = saao_data[:,0] - 2400000.0
times['ECAM'] = ecam_data[:,0] - 2400000.0

fluxes['TESS13'] = tess13_data[:,1]
fluxes['TESS27'] = tess27_data[:,1]
#fluxes['NGTS'] = ngts_data[:,1]
fluxes['SAAO'] = saao_data[:,1]
fluxes['ECAM'] = ecam_data[:,1]

fluxes_error['TESS13'] = tess13_data[:,2]
fluxes_error['TESS27'] = tess27_data[:,2]
#fluxes_error['NGTS'] = ngts_data[:,2]
fluxes_error['SAAO'] = saao_data[:,2]
fluxes_error['ECAM'] = ecam_data[:,2]


##################### INSTRUMENTS #####################
    
lc_instruments_TESS13 = np.empty(len(times['TESS13']), dtype=object)
lc_instruments_TESS13[:] = 'TESS13'
lc_instruments_TESS27 = np.empty(len(times['TESS27']), dtype=object)
lc_instruments_TESS27[:] = 'TESS27'
#lc_instruments_NGTS = np.empty(len(times['NGTS']), dtype=object)
#lc_instruments_NGTS[:] = 'NGTS'
lc_instruments_SAAO = np.empty(len(times['SAAO']), dtype=object)
lc_instruments_SAAO[:] = 'SAAO'
lc_instruments_ECAM = np.empty(len(times['ECAM']), dtype=object)
lc_instruments_ECAM[:] = 'ECAM'

##################### CREATE FINAL INPUT FILES ##################### 

times_all = np.concatenate([np.array(times['TESS13']),np.array(times['TESS27']),
                      np.array(times['SAAO']),np.array(times['ECAM'])])
    
fluxes_all = np.concatenate([np.array(fluxes['TESS13']),np.array(fluxes['TESS27']),
                             np.array(fluxes['SAAO']),np.array(fluxes['ECAM'])])
    
fluxes_error_all = np.concatenate([np.array(fluxes_error['TESS13']),np.array(fluxes_error['TESS27']),
                             np.array(fluxes_error['SAAO']),np.array(fluxes_error['ECAM'])])

lc_instruments_all = np.concatenate([lc_instruments_TESS13,lc_instruments_TESS27,
                                     lc_instruments_SAAO,lc_instruments_ECAM])

np.savetxt('input_data/TOI-5542_input_lc_all.txt',np.column_stack([times_all,fluxes_all,fluxes_error_all,lc_instruments_all]),
            fmt=('%.8f','%.8f','%.8f','%s'))

##################### TESS only ##################### 

times_all = np.concatenate([np.array(times['TESS13']),np.array(times['TESS27'])])
    
fluxes_all = np.concatenate([np.array(fluxes['TESS13']),np.array(fluxes['TESS27'])])
    
fluxes_error_all = np.concatenate([np.array(fluxes_error['TESS13']),np.array(fluxes_error['TESS27'])])

lc_instruments_all = np.concatenate([lc_instruments_TESS13,lc_instruments_TESS27])

#np.savetxt('input_data/TIC466206508_input_lc_TESSonly.txt',np.column_stack([times_all,fluxes_all,fluxes_error_all,lc_instruments_all]),
#            fmt=('%.8f','%.8f','%.8f','%s'))


##################### CREATE FINAL INPUT FILES that need to be added to detrending files ##################### 

times_all = np.concatenate([np.array(times['TESS13']),np.array(times['TESS27'])])
    
fluxes_all = np.concatenate([np.array(fluxes['TESS13']),np.array(fluxes['TESS27'])])
    
fluxes_error_all = np.concatenate([np.array(fluxes_error['TESS13']),np.array(fluxes_error['TESS27'])])

lc_instruments_all = np.concatenate([lc_instruments_TESS13,lc_instruments_TESS27])

#np.savetxt('input_data/TIC466206508_input_lc_detrending.txt',np.column_stack([times_all,fluxes_all,fluxes_error_all,lc_instruments_all]),
#            fmt=('%.8f','%.8f','%.8f','%s'))

############## TESS 13 & 27 GP Files ##############

times_extra = np.concatenate([np.array(times['TESS13']),np.array(times['TESS27'])])

lc_instruments_extra = np.concatenate([lc_instruments_TESS13,lc_instruments_TESS27])
    
np.savetxt('input_data/TOI-5542_input_lc_TESS13gp_TESS27gp.txt',np.column_stack([times_extra,lc_instruments_extra]),
            fmt=('%.8f','%s'))


############## ECAM Detrending ##############

times_ECAM = np.array(times['ECAM'])
fluxes_ECAM = np.array(fluxes['ECAM'])
fluxes_error_ECAM = np.array(fluxes_error['ECAM'])
xshift_ECAM = np.array(ecam_data[:,3])
yshift_ECAM = np.array(ecam_data[:,4])
am_ECAM = np.array(ecam_data[:,5])
fwhm_ECAM = np.array(ecam_data[:,6])
sky_ECAM = np.array(ecam_data[:,7])
exptime_ECAM = np.array(ecam_data[:,8])
amsquared_ECAM = am_ECAM**2
skysquared_ECAM = sky_ECAM**2

#np.savetxt('input_data/TOI-5542_input_lc_ECAM_am_am2_sky_sky2.txt',
#           np.column_stack([times_ECAM,fluxes_ECAM,fluxes_error_ECAM,lc_instruments_ECAM,
#                            am_ECAM,amsquared_ECAM,sky_ECAM,skysquared_ECAM]),
#            fmt=('%.8f','%.8f','%.8f','%s','%.3f','%.3f','%.3f','%.3f'))

### Need to add this to final file to work ####

############## SAAO Detrending ##############

times_SAAO = np.array(times['SAAO'])
fluxes_SAAO = np.array(fluxes['SAAO'])
fluxes_error_SAAO = np.array(fluxes_error['SAAO'])
xshift_SAAO = np.array(saao_data[:,3])
yshift_SAAO = np.array(saao_data[:,4])
am_SAAO = np.array(saao_data[:,5])
fwhm_SAAO = np.array(saao_data[:,6])
sky_SAAO = np.array(saao_data[:,7])
exptime_SAAO = np.array(saao_data[:,8])
amsquared_SAAO = am_SAAO**2
skysquared_SAAO = sky_SAAO**2

    
#### NEED TO COMBINE THESE MANUALLY into final flux file ####
#np.savetxt('input_data/TOI-5542_input_lc_SAAO_am_am2_sky_sky2.txt',
#           np.column_stack([times_SAAO,fluxes_SAAO,fluxes_error_SAAO,lc_instruments_SAAO,
#                            am_SAAO,amsquared_SAAO,sky_SAAO,skysquared_SAAO]),
#            fmt=('%.8f','%.8f','%.8f','%s','%.3f','%.3f','%.3f','%.3f'))


### Need to add this to final file to work ####

# ----------------------------------------------------------------------------
# Prepare input Radial Velocity data data
# ----------------------------------------------------------------------------

rvs_coralie = pd.read_csv('original_data/TIC466206508_COR14_DRS-3-8.rdb', delimiter='\t', usecols=['rjd', 'vrad', 'svrad'])
rvs_coralie = rvs_coralie.drop(0)
rvs_coralie = rvs_coralie.rename(columns={'rjd': 'RJD', 'vrad': 'RV', 'svrad': 'RVerr'})
rvs_coralie.RJD = rvs_coralie.RJD.astype('float')
rvs_coralie.RV = rvs_coralie.RV.astype('float') / 1000.0
rvs_coralie.RVerr = rvs_coralie.RVerr.astype('float') / 1000.0
rvs_coralie['INSTRUMENT'] = 'CORALIE'

rvs_harps = pd.read_csv('original_data/TIC466206508_HARPS15_DRS-3-5.rdb', delimiter='\t', usecols=['rjd', 'vrad', 'svrad'])
rvs_harps = rvs_harps.drop(0)
rvs_harps = rvs_harps.rename(columns={'rjd': 'RJD', 'vrad': 'RV', 'svrad': 'RVerr'})
rvs_harps.RJD = rvs_harps.RJD.astype('float')
rvs_harps.RV = rvs_harps.RV.astype('float') / 1000.0
rvs_harps.RVerr = rvs_harps.RVerr.astype('float') / 1000.0
rvs_harps['INSTRUMENT'] = 'HARPS'

##### REMOVE HARPS IN TRANSIT DATA POINT ##### 2459430.64989303
print(rvs_harps)
rvs_harps = rvs_harps.drop(6)
print(rvs_harps)

# Combining data
rvs_full = rvs_harps.append(rvs_coralie).sort_values(by='RJD')
rvs_filename = 'input_data/TOI-5542_rvs.txt'
rvs_full.to_csv(rvs_filename, index=False, header=False, sep='\t')

#ones_harps = np.ones(len(rvs_harps.RJD))
#ones_coralie = np.ones(len(rvs_coralie.RJD))

#np.savetxt('processed_data/TOI-5542_HARPS.dat',np.column_stack([rvs_harps.RJD-50000,rvs_harps.RV,rvs_harps.RVerr,ones_harps,ones_harps,ones_harps]),
#            fmt=('%.8f','%.8f','%.8f','%i','%i','%i'))

#np.savetxt('processed_data/TOI-5542_CORALIE.dat',np.column_stack([rvs_coralie.RJD-50000,rvs_coralie.RV,rvs_coralie.RVerr,ones_coralie,ones_coralie,ones_coralie]),
#            fmt=('%.8f','%.8f','%.8f','%i','%i','%i'))