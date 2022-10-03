#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Discussion plots for TOI-5542 paper
Created on Thu Sep  1 14:05:40 2022
@author: nolangrieves
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times']})
rc('axes', linewidth=2)
import matplotlib.colors as mcol
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
from scipy import stats
#ignore divide by 0:
np.seterr(divide='ignore', invalid='ignore')

#load planets from NASA Exoplanet Archive
df = pd.read_csv('PS_2022.09.01_05.11.39.csv',sep=',',header=290,low_memory=False)
pnames = np.array(df['pl_name'])

#list(df.columns)

mass = np.array(df['pl_bmassj'])
mass_ue = np.array(abs(df['pl_bmassjerr1']))
mass_le = np.array(abs(df['pl_bmassjerr2']))
masserr = np.array([mass_le,mass_ue]).max(axis=0)

rad = np.array(df['pl_radj'])
rad_ue = np.array(abs(df['pl_radjerr1']))
rad_le = np.array(abs(df['pl_radjerr2']))
raderr = np.array([rad_le,rad_ue]).max(axis=0)

ecc =  np.array(df['pl_orbeccen'])
ecc_ue = np.array(abs(df['pl_orbeccenerr1']))
ecc_le = np.array(abs(df['pl_orbeccenerr1']))
eccerr = np.array([ecc_le,ecc_ue]).max(axis=0)

period = np.array(df['pl_orbper'])
period_ue = np.array(abs(df['pl_orbpererr1']))
period_le = np.array(abs(df['pl_orbpererr2']))
perioderr = np.array([period_le,period_ue]).max(axis=0)

feh = np.array(df['st_met'])
feh_ue = np.array(abs(df['st_meterr1']))
feh_le = np.array(abs(df['st_meterr2']))

dens = np.array(df['pl_dens'])
dens_ue = np.array(abs(df['pl_denserr1']))
dens_le = np.array(abs(df['pl_denserr2']))

strad = np.array(df['st_rad'])
strad_ue = np.array(abs(df['st_raderr1']))
strad_le = np.array(abs(df['st_raderr2']))

stmas = np.array(df['st_mass'])
stmas_ue = np.array(abs(df['st_masserr1']))
stmas_le = np.array(abs(df['st_masserr2']))

stmet = np.array(df['st_met'])
stmet_ue = np.array(abs(df['st_meterr1']))
stmet_le = np.array(abs(df['st_meterr2']))
stmeterr = np.array([stmet_le,stmet_ue]).max(axis=0)

stage = np.array(df['st_age'])
stage_ue = np.array(abs(df['st_ageerr1']))
stage_le = np.array(abs(df['st_ageerr2']))
stageerr = np.array([stage_le,stage_ue]).max(axis=0)

stteff = np.array(df['st_teff'])
stteff_ue = np.array(abs(df['st_tefferr1']))
stteff_le = np.array(abs(df['st_tefferr2']))
sttefferr = np.array([stage_le,stage_ue]).max(axis=0)

vmags = df['sy_vmag']
teff = df['st_teff']
equiltemp = df['pl_eqt']
insol = df['pl_insol']
strad = df['st_rad']
asemi = df['pl_orbsmax']
insol = df['pl_insol']
default_flag = df['default_flag']
disc_ref = df['disc_refname']



################ CREATE GIANT PLANET SAMPLE ################ 
massprecision = 0.25
radiusprecision = 0.08
metalerrorcut = 0.25

####### For periods less the 5 days #######
whg_noecc = ((default_flag == 1)
       & (mass > 0.5) & (mass < 13.0) 
       & (raderr/rad < radiusprecision)
       & (masserr/mass < massprecision) 
       & (stmeterr < metalerrorcut)
       & (np.isfinite(period) == True) 
       & ((period) < 5 ) 
       & (np.isfinite(ecc) == False)
       )

print('Giant planets in sample without eccentricity (P<5 days): ',len(mass[whg_noecc]))

#### ADD in 0 eccentricity
ecc[whg_noecc] = 0
eccerr[whg_noecc] = 0

####### For periods greater than 5 days #######
whg_noecc = ((default_flag == 1)
       & (mass > 0.5) & (mass < 13.0) 
       & (raderr/rad < radiusprecision)
       & (masserr/mass < massprecision) 
       & (stmeterr < metalerrorcut)
       & (np.isfinite(period) == True) 
       & ((period) > 5 ) 
       & (np.isfinite(ecc) == False)
       )

print('Giant planets in sample without eccentricity (P>5 days): ',len(mass[whg_noecc]))

print('-- Giant planets with eccentricity and periods > 5 days --')
print('Planet  Period Ecc Ecc_err')
print('-- -- -- -- -- -- -- -- --')
for i in range(len(pnames[whg_noecc])):
    print(pnames[whg_noecc][i],period[whg_noecc][i],ecc[whg_noecc][i],eccerr[whg_noecc][i])
print('--------------------------')
#Giant planet sample without eccentricity and period >5 days:
#Kepler-56 c 21.40239		ecc = 0.00±0.01
#Kepler-91 b 6.24658		ecc = 0.05±0.02
#WASP-161 b 5.4060425	ecc = 0±0.43
#### Fix individual planet eccentricity values ###

ecc[pnames == 'Kepler-56 c'] = 0
eccerr[pnames == 'Kepler-56 c'] = 0.01
ecc[pnames == 'Kepler-91 b'] = 0.05
eccerr[pnames == 'Kepler-91 b'] = 0.02
ecc[pnames == 'WASP-161 b'] = 0.0
eccerr[pnames == 'WASP-161 b'] = 0.43



whg = ((default_flag == 1)
       & (mass > 0.5) & (mass < 13.0) 
       & (raderr/rad < radiusprecision)
       & (masserr/mass < massprecision) 
       & (stmeterr < metalerrorcut)
       & (np.isfinite(period) == True) 
       & (np.isfinite(ecc) == True)
       )

print('Final Giant planets in sample: ',len(mass[whg]))

a = pnames[whg]
import collections
repeat = [item for item, count in collections.Counter(a).items() if count > 1]

#print(repeat)
print('# REPEATED:',len(repeat))

seen = set()
uniq = []
for x in a:
    if x not in seen:
        uniq.append(x)
        seen.add(x)
        
#print('----------')
#print(seen)
print('# WITH NO REPEATS',len(seen))

############## Add in TIC466206508 / TOI-5542 ##############

period_tic466 = 75.12375
ecc_tic466 = 0.018
eccerr_tic466 = 0.026
mass_tic466 = 1.32
masserr_tic466 = 0.10
rad_tic466 = 1.009
raderr_tic466 = 0.036
stmet_tic466 = -0.21
stmeterr_tic466 = 0.08
stage_tic466 = 10.8
stageerr_tic466 = 3.6
stteff_tic466 = 5700
sttefferr_tic466 = 80
#print(len(mass))

equiltemp_tic466 = 441
insol_tic466 = 9.6

mass = np.append(mass,mass_tic466)
masserr = np.append(masserr,masserr_tic466)
rad = np.append(rad,rad_tic466)
raderr = np.append(raderr,raderr_tic466)
ecc = np.append(ecc,ecc_tic466)
eccerr = np.append(eccerr,eccerr_tic466)
period = np.append(period,period_tic466)
stmet = np.append(stmet,stmet_tic466)
stmeterr = np.append(stmeterr,stmeterr_tic466)

stage = np.append(stage,stage_tic466)
stageerr = np.append(stageerr,stageerr_tic466)

stteff = np.append(stteff,stteff_tic466)
sttefferr = np.append(sttefferr,sttefferr_tic466)

equiltemp = np.append(equiltemp,equiltemp_tic466)
insol = np.append(insol,insol_tic466)

pnames = np.append(pnames,'TOI-5542b')
default_flag = np.append(default_flag,1)

disc_ref = np.append(disc_ref,' target=ref> This Work </a> ')

######### Set new call to grab TIC466206508 ##########

whg = ((default_flag == 1)
       & (mass > 0.5) & (mass < 13.0) 
       & (raderr/rad < radiusprecision)
       & (masserr/mass < massprecision) 
       & (stmeterr < metalerrorcut)
       & (np.isfinite(period) == True) 
       & (np.isfinite(ecc) == True)
       )

print('--------------------------')
print('Final Giant planets in sample with TOI-5542b: ',len(mass[whg]))


##### Double check that doesn't change #####
whg = ((default_flag == 1)
       & (mass > 0.5) & (mass < 13.0) 
       & (raderr/rad < radiusprecision)
       & (masserr/mass < massprecision) 
       & (stmeterr < metalerrorcut)
       & (np.isfinite(period) == True) 
       & (np.isfinite(ecc) == True)
       & (np.isfinite(stteff) == True)
       #& (stageerr/stage < stageprecision) 
       #& (np.isfinite(insol) == True)
       )

print('Final Giant planets in sample with Teff: ',len(mass[whg]))

fout=open('GiantPlanetSample.tex', 'w')
fout.write(r'\begin{table*}'+ "\n")
fout.write(r''+"\n")

# for i in range(len(pnames[whg])):
#     fout.write(str(i+1)+r' & %s & %.2f$\pm$%.2f & %.3f$\pm$%.3f & %.5f & %.3f$\pm$%.3f & .2f$\pm$%.2f & %i & '%(
#         pnames[whg][i],mass[whg][i],masserr[whg][i],rad[whg][i],raderr[whg][i],period[whg][i],
#         ecc[whg][i],eccerr[whg][i],stmet[whg][i],stmeterr[whg][i],stteff[whg][i]))#+str(stage[whg][i])+'$\pm$'+str(stageerr[whg][i])+' & '+disc_ref[whg][i].split('target=ref>')[1].split('</a>')[0]+' \\\ \n')
# fout.write(r''+"\n")
# fout.write(r'\end{table*}'+ "\n")
# fout.close()

stmet = np.array(stmet)
for i in range(len(pnames[whg])):
    fout.write(str(i+1)+r' & %s & %.2f$\pm$%.2f & %.3f$\pm$%.3f & %.5f & %.3f$\pm$%.3f & %.2f$\pm$%.2f & %i & '%(
        pnames[whg][i],mass[whg][i],masserr[whg][i],rad[whg][i],raderr[whg][i],
        period[whg][i],ecc[whg][i],eccerr[whg][i],stmet[whg][i],stmeterr[whg][i],
        stteff[whg][i])+str(stage[whg][i])+'$\pm$'+str(stageerr[whg][i])+' & '+disc_ref[whg][i].split('target=ref>')[1].split('</a>')[0]+' \\\ \n')
fout.write(r''+"\n")
fout.write(r'\end{table*}'+ "\n")
fout.close()


#discovery[i].split('target=ref>')[1].split('</a>')[0]

#################################################
################ Period vs Mass #################
#################################################

fsize = 30
fsize2 = 28

fig, ax = plt.subplots(figsize=(10,8))

plt.scatter(period_tic466, mass_tic466, s=300, facecolors='none', edgecolors='black',linewidth=2)
plt.text(period_tic466+11, mass_tic466,'TOI-5542b',fontsize=20)

img1 = ax.scatter(period[whg],mass[whg],s=120,edgecolors='black',
            c=ecc[whg],cmap='jet',**{"zorder":50},alpha=0.6)#,norm=matplotlib.colors.LogNorm())

cb1 = fig.colorbar(img1,ax=ax)
cb1.set_label(r'Eccentricity',fontsize=fsize)
cb1.ax.tick_params(labelsize=fsize2)

plt.xlabel(r'Period (days)',fontsize=fsize)
plt.ylabel(r'Mass (M$_{\rm{Jup}}$)',fontsize=fsize)
plt.xscale('log')
plt.yscale('log')
plt.xticks([1,10,100,1000],['1','10','100','1000'],fontsize=fsize2)
plt.yticks([1,10],['1','10'],fontsize=fsize2)
ax.tick_params(which='both', width=1.5)
ax.tick_params(which='major', length=12)
ax.tick_params(which='minor', length=8)#, color='r')
ax.tick_params(which='both', direction='in')
ax.tick_params

plt.tight_layout()
plt.savefig('plots/mass_period_ecc.pdf')

#################################################
#############  Rp vs Period with Teff ###########
#################################################

fig, ax = plt.subplots(figsize=(10,8))

plt.scatter(period_tic466, rad_tic466, s=300, facecolors='none', edgecolors='black',linewidth=2)
plt.text(period_tic466, rad_tic466-0.08,'TOI-5542b',fontsize=20)

img1 = ax.scatter(period[whg],rad[whg],s=120,edgecolors='black',
            c=stteff[whg],cmap='jet',**{"zorder":50},alpha=0.6)#,norm=matplotlib.colors.LogNorm())

cb1 = fig.colorbar(img1,ax=ax)
cb1.set_label(r'Stellar T$_{\rm{eff}}$ (K)',fontsize=fsize)
cb1.ax.tick_params(labelsize=fsize2)

plt.xlabel(r'Period (days)',fontsize=fsize)
plt.ylabel(r'Radius (R$_{\rm{Jup}}$)',fontsize=fsize)
plt.xscale('log')
plt.xticks([1,10,100,1000],['1','10','100','1000'],fontsize=fsize2)
plt.yticks(fontsize=fsize2)
ax.tick_params(which='both', width=1.5)
ax.tick_params(which='major', length=12)
ax.tick_params(which='minor', length=8)#, color='r')
ax.tick_params(which='both', direction='in')
ax.tick_params

plt.tight_layout()
plt.savefig('plots/radius_period_teff_all.pdf')


#################################################
#############  Mass vs Metallicity ##############
#################################################

whg_hot = ((default_flag == 1)
           & (mass > 0.5) & (mass < 13.0) 
       & (raderr/rad < radiusprecision)
       & (masserr/mass < massprecision) 
       & (stmeterr < metalerrorcut)
       & (np.isfinite(period) == True) 
       & (np.isfinite(ecc) == True)
       & (period < 10)
       )


whg_warm = ((default_flag == 1)
            & (mass > 0.5) & (mass < 13.0) 
       & (raderr/rad < radiusprecision)
       & (masserr/mass < massprecision) 
       & (stmeterr < metalerrorcut)
       & (np.isfinite(period) == True) 
       & (np.isfinite(ecc) == True)
       & (period >= 10)
       )


#print('Giants with good [Fe/H]: ',len(mass[whg]),' planets')
print('Giant sample with P < 10 (%i planets), mean [Fe/H] = %.5f +/- %.5f'%(len(mass[whg_hot]),np.mean(stmet[whg_hot]),np.std(stmet[whg_hot])/np.sqrt(len(stmet[whg_hot]))))
print('Giant sample with P >= 10 (%i planets), mean [Fe/H] = %.5f +/- %.5f'%(len(mass[whg_warm]),np.mean(stmet[whg_warm]),np.std(stmet[whg_warm])/np.sqrt(len(stmet[whg_warm]))))

fsize=32
fsize2=30

fig = plt.figure(figsize=(12,10))
grid = gridspec.GridSpec(20, 20)
grid.update(wspace=0.0, hspace=0.0) # set the spacing between axes. 
ax = fig.add_subplot(grid[6:,:])
#hist1 = fig.add_subplot(grid[:3,:17])
hist1 = fig.add_subplot(grid[:6,:])


hist1.spines["top"].set_visible(False)  
hist1.spines["right"].set_visible(False)
hist1.spines["left"].set_visible(False)
hist1.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
hist1.xaxis.set_minor_formatter(ticker.NullFormatter())
hist1.set_xticks([])
hist1.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
hist1.yaxis.set_minor_formatter(ticker.NullFormatter())
hist1.set_yticks([])

ax.scatter(stmet[whg_hot],mass[whg_hot],s=120,edgecolors=None,alpha=0.6,color='firebrick',label='P < 10 d')
ax.scatter(stmet[whg_warm],mass[whg_warm],s=120,edgecolors=None,alpha=0.8,color='steelblue',label='P > 10 d')

ax.legend(fontsize=28)

binwidth = 0.1
hmin = -0.7
hmax = 0.7

hist,bin_edges = np.histogram(stmet[whg_hot],bins=np.arange(hmin,hmax,binwidth))
hist = hist#/np.sum(hist)
hist_hot = hist.copy()
hist1.bar((bin_edges[:-1]+binwidth/2.0), hist,width=binwidth,color='firebrick',alpha=0.6,edgecolor='black')#density=True

hist,bin_edges = np.histogram(stmet[whg_warm],bins=np.arange(hmin,hmax,binwidth))
hist = hist#/np.sum(hist)
hist1.bar((bin_edges[:-1]+binwidth/2.0), hist,width=binwidth,color='steelblue',alpha=0.8,edgecolor='black',bottom=hist_hot)

#kde = stats.gaussian_kde(stmet[whg_hot])
#kde_x = np.linspace(hmin, hmax, 1000)
#hist1.plot(kde_x, kde(kde_x),color='firebrick')
#
#kde = stats.gaussian_kde(stmet[whg_warm])
#kde_x = np.linspace(hmin, hmax, 1000)
#hist1.plot(kde_x, kde(kde_x),color='steelblue')

ax.scatter(stmet_tic466, mass_tic466, s=400, facecolors='none', edgecolors='black',linewidth=2)
ax.text(stmet_tic466-0.25, mass_tic466,'TOI-5542b',fontsize=24,color='black')

ax.set_xlim(-0.7,0.7)
hist1.set_xlim(-0.7,0.7)

ax.set_xlabel(r'Stellar Metallicity [Fe/H]',fontsize=fsize)
ax.set_ylabel(r'Mass (M$_{\rm{Jup}}$)',fontsize=fsize)
ax.tick_params(which='both', width=1.5,labelsize=fsize2)
ax.tick_params(which='major', length=12)
ax.tick_params(which='minor', length=8)#, color='r')
ax.tick_params(which='both', direction='in')

ax_t = ax.secondary_xaxis('top')
ax_t.tick_params(axis='x', direction='in')
ax_t.tick_params(axis='x', width=1.5)
ax_t.tick_params(axis='x', length=12)
ax_t.tick_params(axis='x', length=8)

ax.set_yscale('log')
ax.set_yticks([1,10],['1','10'],fontsize=fsize2)

plt.tight_layout()
plt.savefig('plots/mass_metallicity_sep_period.pdf')

#print(pnames[((whg_warm) & (stmet <-0.15))])

whg_lowmass = ((default_flag == 1) &
               (mass > 0.5) & (mass < 4.0) 
       & (raderr/rad < radiusprecision)
       & (masserr/mass < massprecision) 
       & (stmeterr < metalerrorcut)
       & (np.isfinite(period) == True) 
       & (np.isfinite(ecc) == True)
       )

whg_highmass = ((default_flag == 1) &
                (mass >= 4.0) & (mass < 13.0) 
       & (raderr/rad < radiusprecision)
       & (masserr/mass < massprecision) 
       & (stmeterr < metalerrorcut)
       & (np.isfinite(period) == True) 
       & (np.isfinite(ecc) == True)
       )

print('Giant sample: 0.5 Mjup < Mp < 4 Mjup (%i planets) mean [Fe/H] = %.5f +/- %.5f'%(len(stmet[whg_lowmass]),np.mean(stmet[whg_lowmass]),np.std(stmet[whg_lowmass])/np.sqrt(len(stmet[whg_lowmass]))))
print('Giant sample: 4 Mjup < Mp < 13 Mjup (%i planets) mean [Fe/H] = %.5f +/- %.5f'%(len(stmet[whg_highmass]),np.mean(stmet[whg_highmass]),np.std(stmet[whg_highmass])/np.sqrt(len(stmet[whg_highmass]))))


#################################################
#############  Stellar Age ######################
#################################################

stageprecision = 0.4

whg = ((default_flag == 1)
       & (mass > 0.5) & (mass < 13.0) 
       & (raderr/rad < radiusprecision)
       & (masserr/mass < massprecision) 
       & (stmeterr < metalerrorcut)
       & (np.isfinite(period) == True) 
       & (np.isfinite(ecc) == True)
       & (np.isfinite(stage) == True)
       & (stageerr/stage < stageprecision) 
       #& (np.isfinite(insol) == True)
       )

print('Final Giant planets in sample with age: ',len(mass[whg]))


fsize = 30
fsize2 = 28

# fig, ax = plt.subplots(figsize=(10,8))

# plt.scatter(period_tic466, mass_tic466, s=300, facecolors='none', edgecolors='black',linewidth=2)
# plt.text(period_tic466+10, mass_tic466,'TOI-5542b',fontsize=20)

# img1 = ax.scatter(period[whg],mass[whg],s=120,edgecolors='black',
#             c=stage[whg],cmap='jet',**{"zorder":50},alpha=0.6)#,norm=matplotlib.colors.LogNorm())

# cb1 = fig.colorbar(img1,ax=ax)
# cb1.set_label(r'Stellar Age (Gyr)',fontsize=fsize)
# cb1.ax.tick_params(labelsize=fsize2)

# plt.xlabel(r'Period (days)',fontsize=fsize)
# plt.ylabel(r'Mass (M$_{\rm{Jup}}$)',fontsize=fsize)
# plt.xscale('log')
# plt.xticks([1,10,100,1000],['1','10','100','1000'],fontsize=fsize2)
# plt.yticks(fontsize=fsize2)
# ax.tick_params(which='both', width=1.5)
# ax.tick_params(which='major', length=12)
# ax.tick_params(which='minor', length=8)#, color='r')
# ax.tick_params(which='both', direction='in')
# ax.tick_params

# plt.tight_layout()
# #plt.savefig('plots/znot_use/mass_period_age.pdf')

# #################################################

# fig, ax = plt.subplots(figsize=(10,8))

# plt.scatter(period_tic466, rad_tic466, s=300, facecolors='none', edgecolors='black',linewidth=2)
# plt.text(period_tic466, rad_tic466+0.03,'TOI-5542b',fontsize=20)

# img1 = ax.scatter(period[whg],rad[whg],s=120,edgecolors='black',
#             c=stage[whg],cmap='jet',**{"zorder":50},alpha=0.6)#,norm=matplotlib.colors.LogNorm())

# cb1 = fig.colorbar(img1,ax=ax)
# cb1.set_label(r'Stellar Age (Gyr)',fontsize=fsize)
# cb1.ax.tick_params(labelsize=fsize2)

# plt.xlabel(r'Period (days)',fontsize=fsize)
# plt.ylabel(r'Radius (R$_{\rm{Jup}}$)',fontsize=fsize)
# plt.xscale('log')
# plt.xticks([1,10,100,1000],['1','10','100','1000'],fontsize=fsize2)
# plt.yticks(fontsize=fsize2)
# ax.tick_params(which='both', width=1.5)
# ax.tick_params(which='major', length=12)
# ax.tick_params(which='minor', length=8)#, color='r')
# ax.tick_params(which='both', direction='in')
# ax.tick_params

# plt.tight_layout()
# #plt.savefig('plots/znot_use/radius_period_age.pdf')



fig, ax = plt.subplots(figsize=(10,8))

plt.scatter(period_tic466, stage_tic466, s=300, facecolors='none', edgecolors='black',linewidth=2)
plt.text(period_tic466+10, stage_tic466+0.03,'TOI-5542b',fontsize=20)

img1 = ax.scatter(period[whg],stage[whg],s=120,edgecolors='black',
            c=stmet[whg],cmap='jet',**{"zorder":50},alpha=0.6)#,norm=matplotlib.colors.LogNorm())

cb1 = fig.colorbar(img1,ax=ax)
cb1.set_label(r'Stellar Metallicity [Fe/H]',fontsize=fsize)
cb1.ax.tick_params(labelsize=fsize2)

plt.xlabel(r'Period (days)',fontsize=fsize)
plt.ylabel(r'Stellar Age (Gyr)',fontsize=fsize)
plt.xscale('log')
plt.xticks([1,10,100,1000],['1','10','100','1000'],fontsize=fsize2)
plt.yticks(fontsize=fsize2)
ax.set_ylim(0,13.5)
ax.tick_params(which='both', width=1.5)
ax.tick_params(which='major', length=12)
ax.tick_params(which='minor', length=8)#, color='r')
ax.tick_params(which='both', direction='in')
ax.tick_params

plt.tight_layout()
plt.savefig('plots/age_period_metal.pdf')




# #################################################

# #5000-6000
# whg = ((default_flag == 1) 
#        & (mass > 0.5) & (mass < 13.0) 
#        & (raderr/rad < radiusprecision)
#        & (masserr/mass < massprecision) 
#        & (stmeterr < metalerrorcut)
#        & (np.isfinite(period) == True) 
#        & (np.isfinite(ecc) == True)
#        & (stteff > 5000) & (stteff < 6000)
#        #& (stageerr/stage < stageprecision) 
#        #& (np.isfinite(insol) == True)
#        )

# print('Final Giant planets in sample with Teff 5000-6000: ',len(mass[whg]))


# fig, ax = plt.subplots(figsize=(10,8))

# plt.scatter(period_tic466, rad_tic466, s=300, facecolors='none', edgecolors='black',linewidth=2)
# plt.text(period_tic466, rad_tic466-0.08,'TOI-5542b',fontsize=20)

# img1 = ax.scatter(period[whg],rad[whg],s=120,edgecolors='black',
#             c=stteff[whg],cmap='jet',**{"zorder":50},alpha=0.6)#,norm=matplotlib.colors.LogNorm())

# cb1 = fig.colorbar(img1,ax=ax)
# cb1.set_label(r'Stellar T$_{\rm{eff}}$ (K)',fontsize=fsize)
# cb1.ax.tick_params(labelsize=fsize2)

# plt.xlabel(r'Period (days)',fontsize=fsize)
# plt.ylabel(r'Radius (R$_{\rm{Jup}}$)',fontsize=fsize)
# plt.xscale('log')
# plt.xticks([1,10,100,1000],['1','10','100','1000'],fontsize=fsize2)
# plt.yticks(fontsize=fsize2)
# ax.tick_params(which='both', width=1.5)
# ax.tick_params(which='major', length=12)
# ax.tick_params(which='minor', length=8)#, color='r')
# ax.tick_params(which='both', direction='in')
# ax.tick_params

# plt.tight_layout()
# #plt.savefig('plots/znot_use/radius_period_teff_Gdwarfs.pdf')


