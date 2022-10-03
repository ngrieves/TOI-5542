#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 11:14:33 2022

@author: nolangrieves
"""

from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip
import numpy as np
from tqdm import tqdm

###############################################################################
###############################################################################
# Load and plot the data
###############################################################################
###############################################################################
t = Table.read('original_data/TIC-466206508_NGTS_full.fits')
time_offset=int(np.floor(np.min(t['BJD'])))
f = plt.figure(figsize=(15,5))
plt.scatter(t['BJD']-time_offset, t['TARGET_DETRENDED'], c='k', s=1, label = 'Raw NGTS data')
plt.xlabel('BJD - {:,}'.format(time_offset))
plt.ylabel('Detrended Flux')
plt.legend()

###############################################################################
###############################################################################
# OK, lets sigma clip this as there are some serious  outliers
###############################################################################
###############################################################################
mask = (t['TARGET_DETRENDED'] > 0.015) & (t['TARGET_DETRENDED'] < 0.02) & (t['TARGET_DETRENDED_ERR'] < 0.0002)

t = Table.read('original_data/TIC-466206508_NGTS_full.fits')[mask]
f = plt.figure(figsize=(15,5))
plt.scatter(t['BJD'] - time_offset, t['TARGET_DETRENDED'], c='k', s=1, label = 'Masked NGTS data')
plt.xlabel('BJD - {:,}'.format(time_offset))
plt.ylabel('Detrended Flux')
plt.legend()

###############################################################################
###############################################################################
# OK, now lets get the cameras for each night
###############################################################################
###############################################################################

#nights = bruce.data.find_nights_from_data(t['BJD'], dx_lim=0.2)
nights = np.unique(t['DATE'])
#nights = np.unique(np.around(t['BJD']))
print('Found {:,} nights'.format(len(nights)))

unique_cameras = np.unique(np.array(t['CAMERA']))
unique_cameras_count = np.zeros_like(unique_cameras)
unique_cameras_mask = np.zeros((unique_cameras.shape[0], len(nights)), dtype = bool)

camera_per_night = []
for i in tqdm(range(len(nights))):
    #uniques = np.unique(np.array(t[nights[i]]['CAMERA']))
    uniques = np.unique(np.array(t[t['DATE'] == nights[i]]['CAMERA']))
    camera_per_night.append(uniques)
    for j in range(len(uniques)):
        for k in range(len(unique_cameras)):
            if uniques[j]==unique_cameras[k] : 
                unique_cameras_count[k] +=1
                unique_cameras_mask[k, i] = True

f = plt.figure(figsize=(15,5))
for i in range(len(nights)):
    plt.scatter(i*np.ones(len(camera_per_night[i])), camera_per_night[i], c='k')
plt.xlabel('Night')
plt.ylabel('Camera ID')


f = plt.figure()
plt.imshow(unique_cameras_mask, aspect='auto', origin='lower', interpolation=None)
plt.xlabel('Night')
plt.ylabel('Camera ID')

###############################################################################
###############################################################################
# Now lets align each camera's median values based, looking at the best ones that have already been matched first
###############################################################################
###############################################################################
# Lets do a dry run to see if we can find offsets
max_iters = len(unique_cameras)
count = 0
medians = np.zeros(len(unique_cameras))
camera_sort = np.argsort(unique_cameras_count)[::-1]
unique_cameras = unique_cameras[camera_sort]
unique_cameras_count = unique_cameras_count[camera_sort]
unique_cameras_mask = unique_cameras_mask[camera_sort]

count = 0
nnigh_align = 5
while count < len(unique_cameras):
    if count == -1:
        # This is the start
        medians[np.argmax(unique_cameras_count)] = np.median(t['TARGET_DETRENDED'][t['CAMERA']==unique_cameras[np.argmax(unique_cameras_count)]])
    
    cams_without_offset = np.where(medians==0)[0]
    cams_with_offset = np.where(medians>0)[0]

    for i in range(cams_without_offset.shape[0]):
        for j in range(cams_with_offset.shape[0]):
            if np.sum(unique_cameras_mask[cams_without_offset[i]] & unique_cameras_mask[cams_with_offset[j]]) > nnigh_align:
                medians[cams_without_offset[i]] = np.median(t['TARGET_DETRENDED'][t['CAMERA']==unique_cameras[cams_without_offset[i]]])
                medians[cams_without_offset[i]] += (medians[cams_without_offset[i]] - medians[cams_with_offset[j]])
                break
                
    cams_without_offset = np.where(medians==0)[0]
    cams_with_offset = np.where(medians>0)[0]
    if cams_without_offset.shape[0]==0:
        print('All cameras aligned.')
        break
    else:
        print('Cameras {:} havent been aligned yet after iteration {:}'.format(unique_cameras[cams_without_offset], count))
                
    count +=1
    
print('Cam Count        Offset')
for i in range(len(unique_cameras)):
    print('{:} [{:,} nights] : {:.5f}'.format(unique_cameras[i], unique_cameras_count[i], medians[i]))
    
# Lets do a dry run to see if we can find offsets
max_iters = len(unique_cameras)
count = 0
medians = np.zeros(len(unique_cameras))+1
print(medians)
camera_sort = np.argsort(unique_cameras_count)[::-1]
unique_cameras = unique_cameras[camera_sort]
unique_cameras_count = unique_cameras_count[camera_sort]
unique_cameras_mask = unique_cameras_mask[camera_sort]
unique_cams_offset = np.zeros(len(unique_cameras), dtype = bool)
count = 0
nnigh_align = 5
while count < len(unique_cameras):    
    if count ==0 : 
        medians[0] = np.median(t['TARGET_DETRENDED'][t['CAMERA']==unique_cameras[0]])
        count +=1
        unique_cams_offset[0] = True
        continue
        
    cams_without_offset = np.where(~unique_cams_offset)[0]
    cams_with_offset = np.where(unique_cams_offset)[0]

    for i in range(cams_without_offset.shape[0]):
        for j in range(cams_with_offset.shape[0]):
            print('a')


            if np.sum(unique_cameras_mask[cams_without_offset[i]] & unique_cameras_mask[cams_with_offset[j]]) > nnigh_align:
                medians[cams_without_offset[i]] = np.median(t['TARGET_DETRENDED'][t['CAMERA']==unique_cameras[cams_without_offset[i]]])
                medians[cams_without_offset[i]] = (medians[cams_without_offset[i]] - medians[cams_with_offset[j]])
                print('b')
                break
                
    cams_without_offset = np.where(~unique_cams_offset)[0]
    cams_with_offset = np.where(unique_cams_offset)[0]
    if cams_without_offset.shape[0]==0:
        print('All cameras aligned.')
        break
    else:
        print('Cameras {:} havent been aligned yet after iteration {:}'.format(unique_cameras[cams_without_offset], count))
                
    count +=1
    
print('Cam Count        Offset')
for i in range(len(unique_cameras)):
    print('{:} [{:,} nights] : {:.5f}'.format(unique_cameras[i], unique_cameras_count[i], medians[i]))

###############################################################################
###############################################################################
# Now lets detrend
###############################################################################
###############################################################################
t['detrended'] = np.zeros(len(t))
t['detrended_err'] = np.zeros(len(t))
for i in range(len(unique_cameras)):
    mask = t['CAMERA']==unique_cameras[i]
    if i > 0:
        t['detrended'][mask] = t['TARGET_DETRENDED'][mask] / (medians[0] + medians[i])
        t['detrended_err'][mask] = t['TARGET_DETRENDED_ERR'][mask] / (medians[0] + medians[i])
    else: 
        t['detrended'][mask] = t['TARGET_DETRENDED'][mask] / (medians[0])
        t['detrended_err'][mask] = t['TARGET_DETRENDED_ERR'][mask] / (medians[0])        
        
f = plt.figure(figsize=(15,5))
plt.scatter(t['BJD'] - time_offset, t['detrended'], c='k', s=1)
plt.xlabel('BJD - {:}'.format(time_offset))
plt.ylabel('Detrended Flux')

###############################################################################
###############################################################################
# OK, lets bin this and see if anyt offset exists
###############################################################################
###############################################################################
def fancy_bin(time, data, error, binw):
	"""
	Bin input arrays to bwin. Arrays do not have to be sorted.
	Data is averaged werighted by the inverse error.
	"""
	start = np.min(time)
	mid = start+binw/2.
	end = start+binw
	
	bmids = []
	bvals = []
	berror = []
	bin_n = 0	
	
	while start < np.max(time):
		bin_n+=1
		n = 0
		tot = 0
		lims = np.where(np.logical_and(time<end,time>=start))
		n = len(data[lims])
		tot = np.sum(x / y for x, y in zip(data[lims], error[lims])) / np.sum(1.0/error[lims]) #np.sum is depracated
		#tot = np.sum(x / y for x, y in zip(data[lims], error[lims])) / sum(1.0/error[lims]) 
		tot_error = sum(error[lims]**2)
		
		if n > 0:
			bmids.append(mid)
			bvals.append(tot)
			berror.append(np.sqrt(tot_error)/n)
		### Debugging	
#		print("phase : %.2f - %.2f, mid: %.2f... %.2f +/- %.2f. #obs = %s" % (start, end, mid, tot, sqrt(tot_error)/n, n))
#		print(data[lims], error[lims])
		start +=binw
		mid +=binw
		end +=binw
	
	bmids = np.array(bmids)
	bvals = np.array(bvals)
	beror = np.array(berror)
	return(bmids, bvals, berror)

def bin_data(time, flux, bin_width):
    '''
    Function to bin the data into bins of a given width. time and bin_width 
    must have the same units
    '''

    edges = np.arange(np.min(time), np.max(time), bin_width)
    dig = np.digitize(time, edges)
    time_binned = (edges[1:] + edges[:-1]) / 2
    flux_binned = np.array([np.nan if len(flux[dig == i]) == 0 else flux[dig == i].mean() for i in range(1, len(edges))])
    err_binned = np.array([np.nan if len(flux[dig == i]) == 0 else sem(flux[dig == i]) for i in range(1, len(edges))])
    time_bin = time_binned[~np.isnan(err_binned)]
    err_bin = err_binned[~np.isnan(err_binned)]
    flux_bin = flux_binned[~np.isnan(err_binned)]   

    return time_bin, flux_bin, err_bin

#t_bin, f_bin, fe_bin = bruce2.data.bin_data(t['BJD'], t['detrended'], 0.6)
#binwidth = 30.0/60.0/24.0 #in days
binwidth = 30.0/60.0/24.0 #in days
t_bin, f_bin, fe_bin = fancy_bin(t['BJD'], t['detrended'],t['detrended_err'],binwidth)
#fe_bin = np.array(fe_bin)

f = plt.figure(figsize=(15,5))
plt.scatter(t['BJD'] - time_offset, t['detrended'], c='k', s=1, label = 'NGTS raw')
plt.scatter(t_bin - time_offset, f_bin, c='b', s=1, label = 'NGTS 30 minute bin')
plt.xlabel('BJD - {:}'.format(time_offset))
plt.ylabel('Detrended Flux')
plt.legend()

###############################################################################
###############################################################################
# plot the binned data
###############################################################################
###############################################################################

fig, (ax1, ax2) = plt.subplots(2, 1,sharex=False,
     gridspec_kw={'height_ratios': [1, 0.9]},figsize=(15,10))
ax1.errorbar(t_bin-time_offset, f_bin,yerr=fe_bin, fmt='.', color='black',alpha=0.6,linewidth=1, label = 'NGTS 30 minute bin')

T0 = 2459430.58577
ax1.axvline(T0-time_offset,linestyle='--',color='red',alpha=0.5,linewidth=2,label='3 August 2021')
#plt.ylim(0.99, 1.02)
ax1.set_ylim(0.98, 1.02)
ax1.legend(fontsize=22)
ax1.set_xlabel('BJD - {:}'.format(time_offset),fontsize=28)
ax1.set_ylabel('Detrended Flux',fontsize=28)

from astropy.timeseries import LombScargle
from scipy import signal
max_period = 400
min_period = 0.9

min_frequency = 1./max_period
max_frequency = 1./min_period

ls = LombScargle(t_bin,f_bin,fe_bin)
#freq, power = ls.autopower(minimum_frequency=0.0001, maximum_frequency=1)
freq, power = ls.autopower(minimum_frequency=min_frequency, maximum_frequency=max_frequency)
period = 1./freq
print(period[np.argmax(power)])
print(period[np.argsort(power)[-2]])

#FAP
ls.false_alarm_probability(power.max())  
probabilities = [0.1, 0.01, 0.001]
faps = ls.false_alarm_level(probabilities)

print(np.argsort(power)[1])
ax2.plot(period, power,color='cornflowerblue',linewidth=2,label='NGTS periodogram')
#ax2.axhline(y=faps[2], color='firebrick',linewidth=3, linestyle='--',label='0.1% FAP')

ax2.set_xscale('log')
ax2.set_ylabel('Normalized power',fontsize=28)
ax2.set_xlabel('Period (days)',fontsize=28)
ax2.legend(fontsize=22)
ax2.set_xticks([1,10,100])
ax2.set_xticklabels(['1','10','100'])
ax1.tick_params(axis='both', which='major', labelsize=24,width=3,length=10)
ax2.tick_params(axis='both', which='major', labelsize=24,width=3,length=10)
ax2.tick_params(axis='both', which='minor',width=1.5,length=5)

plt.tight_layout()
plt.savefig('plots/TOI-5542_NGTS_30min.pdf')