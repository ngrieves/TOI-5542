#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Radial Velocity & transit model for TOI-5542 with Juliet
Created on Thu Sep  1 11:26:03 2022
@author: nolangrieves
"""

import juliet
import matplotlib.pyplot as plt
import numpy as np
import corner
import matplotlib.gridspec as gridspec
import pykima as pk
from juliet_definitions import *
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times']})\


star = 'TOI-5542'
output_folder = 'juliet_output_TOI-5542_1sep2022'

input_lc_file = 'data/input_data/TOI-5542_input_lc_all.txt'
input_rv_file = 'data/input_data/TOI-5542_rvs.txt'
input_GP_lc_file = 'data/input_data/TOI-5542_input_lc_TESS13gp_TESS27gp.txt'

####### Stellar Parameters #######
Mstar = 0.890 #M_sun
Mstar_err = 0.056 #M_sun
Rstar = 1.058 #R_sun
Rstar_err = 0.036 #R_sun

teff = 5700
logg = 4.2

rho_star = 1.07 * 1000.0 
rho_star_err = 0.13 * 1000.0 

####### Initial Parameters #######
Tc1_initial = 58679.3
P1_initial = 75.12

# ----------------------------------------------------------------------------
# Priors
# ----------------------------------------------------------------------------
# Define the photometry priors
PRIORS = [
    ('P_p1', 'uniform',[P1_initial-0.2,P1_initial+0.2]),# ('P_p1', 'uniform', [10.0,30.0]),
    ('t0_p1', 'uniform',[Tc1_initial-0.2, Tc1_initial+0.2]), #'uniform', [Tc1_initial-0.1,Tc1_initial+0.1]),
    #('p_p1', 'uniform', [0.0,1.0]),
    #('b_p1', 'uniform', [0.0,1.0]), #('a_p1', 'loguniform', [1,1000.0]),
    ('r1_p1', 'uniform', [0.0,1.0]),
    ('r2_p1', 'uniform', [0.0,1.0]), #('a_p1', 'loguniform', [1,1000.0]),    
    ('q1_TESS13_TESS27', 'normal',[0.335461,0.011]), #from Adrien Deline modification of Espinoza code
    ('q2_TESS13_TESS27', 'normal',[0.260712,0.030]), #from Adrien Deline modification of Espinoza code
    #('q1_NGTS', 'fixed',0.390872), #from Adrien Deline modification of Espinoza code
    #('q2_NGTS', 'fixed',0.281733), #from Adrien Deline modification of Espinoza code
    ('q1_SAAO', 'fixed',0.520824), #from Adrien Deline modification of Espinoza code
    ('q2_SAAO', 'fixed',0.324503), #from Adrien Deline modification of Espinoza code
    ('q1_ECAM', 'fixed',0.520824), #from Adrien Deline modification of Espinoza code
    ('q2_ECAM', 'fixed',0.324503), #from Adrien Deline modification of Espinoza code
    #('ecc_p1', 'uniform', [0.0,1.0]), 
    #('ecc_p1', 'uniform', [0.0,0.4]), 
    ('ecc_p1', 'beta', [0.867,3.03]),
    ('omega_p1', 'uniform', [0.0,360.0]),
    #('esinomega_p1', 'uniform', [-1.0,1.0]),
    #('ecosomega_p1', 'uniform', [-1.0,1.0]),
    ('rho', 'normal',[rho_star,rho_star_err]),
    ('mdilution_TESS13', 'fixed', 1),
    ('mdilution_TESS27', 'fixed', 1),
    #('mdilution_NGTS_SAAO_ECAM', 'fixed', 1.0),
    ('mdilution_SAAO_ECAM', 'fixed', 1.0),
    ('mflux_TESS13', 'normal', [0.0,0.01]),
    #('sigma_w_TESS13', 'loguniform', [1e-5,100.0]),
    ('sigma_w_TESS13', 'fixed', 0),
    ('mflux_TESS27', 'normal', [0.0,0.01]),
    #('sigma_w_TESS27', 'loguniform', [1e-5,100.0]),
    ('sigma_w_TESS27', 'fixed', 0),
    #('mflux_NGTS', 'normal', [0.0,0.01]),
    #('sigma_w_NGTS', 'loguniform', [1e-5,10000.0]),
    ('mflux_SAAO', 'normal', [0.0,0.01]),
    ('sigma_w_SAAO', 'loguniform', [1e-5,10000.0]),
    ('mflux_ECAM', 'normal', [0.0,0.01]),
    ('sigma_w_ECAM', 'loguniform', [1e-5,10000.0]),
    #('theta0_ECAM', 'uniform', [-1,1]),
    #('theta1_ECAM', 'uniform', [-1,1]),
    #('theta2_ECAM', 'uniform', [-1,1]),
    #('theta3_ECAM', 'uniform', [-1,1]),
    #('theta0_SAAO', 'uniform', [-1,1]),
    #('theta1_SAAO', 'uniform', [-1,1]),
    #('theta2_SAAO', 'uniform', [-1,1]),
    #('theta3_SAAO', 'uniform', [-1,1]),
    #('GP_sigma_NGTS', 'loguniform', [1e-6, 1]),
    #('GP_rho_NGTS', 'loguniform', [1e-6,1e3]),
    ('GP_sigma_TESS13', 'loguniform', [1e-6, 1]),
    ('GP_rho_TESS13', 'loguniform', [1e-6,1e3]),
    ('GP_sigma_TESS27', 'loguniform', [1e-6, 1]),
    ('GP_rho_TESS27', 'loguniform', [1e-6,1e3])
]

priors = {}

for prior in PRIORS:
    param = prior[0]
    priors[param] = {}
    priors[param]['distribution'], priors[param]['hyperparameters'] = prior[1], prior[2]   


#-----------------------------------------------------------------------------------
# Define the RV priors
# Radial velocitie in km/s
PRIORSRV = [
    ('K_p1', 'uniform',[0.,1.]),
    ('mu_CORALIE', 'uniform', [-100, 100]),
    ('sigma_w_CORALIE', 'loguniform', [1e-5, 0.2]),
    ('mu_HARPS', 'uniform', [-100, 100]),
    ('sigma_w_HARPS', 'loguniform', [1e-5, 0.2]),
]

# Populate the priors dictionary:
for prior in PRIORSRV:
    param = prior[0]
    priors[param] = {}
    priors[param]['distribution'], priors[param]['hyperparameters'] = prior[1], prior[2]
    
#################################################################
########################## FIT DATA #############################
#################################################################
# Loading the dataset
dataset = juliet.load(priors=priors, 
                      lcfilename = input_lc_file,
                      rvfilename=input_rv_file,
                      GPlceparamfile=input_GP_lc_file,
                      out_folder=output_folder,
                      ld_laws='quadratic')

# Fit dataset
results = dataset.fit(n_live_points=1000,sampler='dynesty',delta_z_lim=0.1) #,nthreads=4)


#################################################################
######################## ANALYZE THE FIT ########################
#################################################################

########################################################
################ PLOT Posterior CORNERS ################
########################################################

# CORNER planetary parameters
first_time = True
posterior_names = ['P_p1', 't0_p1','r1_p1','r2_p1', 'ecc_p1', 'K_p1', 'rho', 'omega_p1','q1_TESS13_TESS27','q2_TESS13_TESS27']
for param_name in posterior_names:
    print(param_name)
    if first_time:
        posterior_data = results.posteriors['posterior_samples'][param_name]
        first_time = False
    else:
        posterior_data  = np.vstack((posterior_data,results.posteriors['posterior_samples'][param_name]))
posterior_data = posterior_data.T

figure = corner.corner(posterior_data, labels=posterior_names)
plt.tight_layout()
plt.savefig(output_folder+'/corner_planet_'+star+'.pdf')
figure.show()

# CORNER intrument
first_time = True
posterior_names = ['mu_HARPS', 'sigma_w_HARPS',
                   'mu_CORALIE', 'sigma_w_CORALIE',
                   'mflux_TESS13',# 'sigma_w_TESS13',
                   'mflux_TESS27',# 'sigma_w_TESS27',
                   #'mflux_NGTS', 'sigma_w_NGTS',
                   'mflux_SAAO', 'sigma_w_SAAO',
                   'mflux_ECAM', 'sigma_w_ECAM',
                   #'GP_sigma_NGTS','GP_rho_NGTS',
                   'GP_sigma_TESS13','GP_rho_TESS13',
                   'GP_sigma_TESS27','GP_rho_TESS27']#,
                   #'theta0_ECAM','theta1_ECAM','theta2_ECAM','theta3_ECAM',
                   #'theta0_SAAO','theta1_SAAO','theta2_SAAO','theta3_SAAO']
first_time = True
for param_name in posterior_names:
    print(param_name)
    if first_time:
        posterior_data = results.posteriors['posterior_samples'][param_name]
        first_time = False
    else:
        posterior_data  = np.vstack((posterior_data,results.posteriors['posterior_samples'][param_name]))
posterior_data = posterior_data.T

figure = corner.corner(posterior_data, labels=posterior_names)
plt.tight_layout()
plt.savefig(output_folder+'/corner_instrument_'+star+'.pdf')
figure.show()


########################################################
################ Get Parameters from posteriors ########
########################################################


G = 6.67408e-11
try:
    a = ((results.posteriors['posterior_samples']['rho']
          *G*((results.posteriors['posterior_samples']['P_'+'p1']*24.*3600.)**2))/(3.*np.pi))**(1./3.)
except KeyError:
    a = results.posteriors['posterior_samples']['a_p1']

# fitted params
P = results.posteriors['posterior_samples']['P_p1']  # orbital period
t0 = results.posteriors['posterior_samples']['t0_p1']

r1 = results.posteriors['posterior_samples']['r1_p1']  
r2 = results.posteriors['posterior_samples']['r2_p1']  
b,p = juliet.utils.reverse_bp(r1, r2, 0., 1.) # impact parameter b = (a / Rs) *cos(i) # planet to star ratio k = (Rp / Rs)
#b = results.posteriors['posterior_samples']['b_p1']  # impact parameter b = (a / Rs) *cos(i)
#p = results.posteriors['posterior_samples']['p_p1']  # planet to star ratio k = (Rp / Rs)

# get inclination from eccentricity and omega from Juliet
ecc = results.posteriors['posterior_samples']['ecc_p1']
omega = results.posteriors['posterior_samples']['omega_p1']

ecc_factor = (1. + ecc*np.sin(omega))/(1. - ecc**2)
inc_inv_factor = (b/a)*ecc_factor
inc = np.arccos(inc_inv_factor)*180./np.pi

# Time of transit to time of periastron
tp = timetrans_to_timeperi(t0, per=p, ecc=ecc, omega=omega*np.pi/180.0)

Rstar_array = np.random.normal(Rstar,Rstar_err,len(p)) * 9.73116 #Rsun to Rjupiter 109.076 #Rsun to Rearth
radius_array = p*Rstar_array  # planet to star ratio k = (Rp / Rs)

radius_value = (get_quantiles(radius_array)[0].round(3))
radius_value_up = round((get_quantiles(radius_array)[1]) - (get_quantiles(radius_array)[0]), 3)
radius_value_low = abs(round((get_quantiles(radius_array)[2]) - (get_quantiles(radius_array)[0]), 3))
t_rad = r"Rp = $%s^{+%s}_{-%s} R_{\oplus}$" % (radius_value, radius_value_up, radius_value_low)

t0_value = get_quantiles(t0)[0].round(4)
t0_value_up = round(get_quantiles(t0)[1] - get_quantiles(t0)[0], 4)
t0_value_low = abs(round(get_quantiles(t0)[2] - get_quantiles(t0)[0], 4))
t_t0 = r'T0 = $%s^{+%s}_{-%s}$' % (t0_value, t0_value_up, t0_value_low)

depth_value = ((get_quantiles(p)[0])**2).round(4)  # depth = p **2 = (Rp/Rs)**2
depth_value_up = round(((get_quantiles(p)[1])**2) - ((get_quantiles(p)[0])**2), 4)
depth_value_low = abs(round(((get_quantiles(p)[2])**2) - ((get_quantiles(p)[0])**2), 4))
t_depth = r'depth = $%s^{+%s}_{-%s}$' % (depth_value, depth_value_up, depth_value_low)

b_value = (get_quantiles(b)[0].round(2))
b_value_up = round((get_quantiles(b)[1]) - (get_quantiles(b)[0]), 2)
b_value_low = abs(round((get_quantiles(b)[2]) - (get_quantiles(b)[0]), 2))
t_b = r'b = $%s^{+%s}_{-%s}$' % (b_value, b_value_up, b_value_low)
                  
ecc_value = (get_quantiles(ecc)[0].round(2))
ecc_value_up = round((get_quantiles(ecc)[1]) - (get_quantiles(ecc)[0]), 2)
ecc_value_low = abs(round((get_quantiles(ecc)[2]) - (get_quantiles(ecc)[0]), 2))
t_ecc = r'b = $%s^{+%s}_{-%s}$' % (ecc_value, ecc_value_up, ecc_value_low)

ecc_value = (get_quantiles(ecc)[0].round(2))
ecc_value_up = round((get_quantiles(ecc)[1]) - (get_quantiles(ecc)[0]), 2)
ecc_value_low = abs(round((get_quantiles(ecc)[2]) - (get_quantiles(ecc)[0]), 2))
t_ecc = r'ecc = $%s^{+%s}_{-%s}$' % (ecc_value, ecc_value_up, ecc_value_low)

per_value = (get_quantiles(P)[0].round(6))
per_value_up = round((get_quantiles(P)[1]) - (get_quantiles(P)[0]), 6)
per_value_low = abs(round((get_quantiles(P)[2]) - (get_quantiles(P)[0]), 6))
t_per = r'P = $%s^{+%s}_{-%s}$ days' % (per_value, per_value_up, per_value_low)

inc_value = (get_quantiles(inc)[0].round(2))
inc_value_up = round((get_quantiles(inc)[1]) - (get_quantiles(inc)[0]), 2)
inc_value_low = abs(round((get_quantiles(inc)[2]) - (get_quantiles(inc)[0]), 2))
t_inc = r'inc = $%s^{+%s}_{-%s}$ deg' % (inc_value, inc_value_up, inc_value_low)

# Planet MASS
Mpsini = pk.utils.get_planet_mass(results.posteriors['posterior_samples']['P_p1'],
                                  results.posteriors['posterior_samples']['K_p1']*1000.0,
                                  e=results.posteriors['posterior_samples']['ecc_p1'],
                                  star_mass=np.random.normal(Mstar,Mstar_err,len(results.posteriors['posterior_samples']['P_p1'])), full_output=True)

Mp = Mpsini[2] / np.sin((inc*np.pi/180.))
Mp_value = round(get_quantiles(Mp)[0], 2)
Mp_value_up = round((get_quantiles(Mp)[1]) - (get_quantiles(Mp)[0]), 2)
Mp_value_low = abs(round((get_quantiles(Mp)[2]) - (get_quantiles(Mp)[0]), 2))
t_Mp = r"Mp = $%s^{+%s}_{-%s} M_{\oplus}$" % (Mp_value, Mp_value_up, Mp_value_low)

K = results.posteriors['posterior_samples']['K_p1'] * 1000   # m/s
K_value = (get_quantiles(K)[0].round(2))
K_value_up = round((get_quantiles(K)[1]) - (get_quantiles(K)[0]), 2)
K_value_low = abs(round((get_quantiles(K)[2]) - (get_quantiles(K)[0]), 2))
t_K = r'K = $%s^{+%s}_{-%s}$ m/s' % (K_value, K_value_up, K_value_low)

#stellar density in kg/m3
rho = results.posteriors['posterior_samples']['rho']
rho_value = get_quantiles(rho)[0].round(4)
rho_value_up = round(get_quantiles(rho)[1] - get_quantiles(rho)[0], 4)
rho_value_low = abs(round(get_quantiles(rho)[2] - get_quantiles(rho)[0], 4))
#change to g/cm3
rho_value = rho_value/1000.0
logg_rho = np.log10(rho_value)


# TRANSIT DURATION
# Transit durations approximations (eq. 14, 15, 16 from Winn 2014)
#k = results.posteriors['posterior_samples']['K_p1'] #* 1000   # m/s
#ecc = results.posteriors['posterior_samples']['ecc_p1']
#omega = results.posteriors['posterior_samples']['omega_p1']
#b,p = juliet.utils.reverse_bp(r1, r2, 0., 1.)
#P = results.posteriors['posterior_samples']['P_p1']  # orbital period

ec_factor = np.sqrt((1. - ecc*ecc)) / (1.0 + ecc*np.sin(omega*np.pi/180.0))

#total transit duration (from contact 1 to 4)
trt = np.sqrt((1. + K/1000)**2 - b** 2) / (a * np.sin(inc*np.pi/180.0))
trt = P / np.pi * np.arcsin(trt) * ec_factor * 24.0
print('total transit duration (from contact 1 to 4): ', round(np.nanmean(trt), 4), 'hours')

#full transit duration (from contact 2 to 3)
tri = np.sqrt((1. - K/1000)**2 - b** 2) / (a * np.sin(inc*np.pi/180))
tri = P / np.pi * np.arcsin(tri) * ec_factor * 24.0
print('full transit duration (from contact 2 to 3): ', round(np.nanmean(tri), 4), 'hours')

#t = 'Juliet params:\n\n' + t_per + '\n' + t_rad + '\n' + t_Mp + '\n' + t_t0 + '\n' + t_depth + '\n' + t_dur + '\n' + t_tau + '\n' + t_b + '\n' + t_ecc + '\n' + t_inc + '\n' + t_K
t = 'Juliet params:\n\n' + t_per + '\n' + t_rad + '\n' + t_Mp + '\n' + t_t0 + '\n' + t_depth + '\n' + t_b + '\n' + t_ecc + '\n' + t_inc + '\n' + t_K
plt.figure()
plt.text(0.1, 0.5, t)
plt.xticks([])
plt.yticks([])
plt.axis('off')
plt.tight_layout()
plt.savefig(output_folder+'/parameters_'+star+'.pdf')



#################################################################
############# PLOT FULL RVS and residuals  ######################
#################################################################
rc('axes', linewidth=4)

# Plot RV times
min_time, max_time = np.min(dataset.times_rv['CORALIE'])-5, np.max(dataset.times_rv['CORALIE'])+5
model_rv_times = np.linspace(min_time, max_time, 10000)

# Evaluate RV model use all the posterior samples also extract model components
rv_model, median_model, components = results.rv.evaluate('HARPS', t=model_rv_times,
                                                         nsamples=150,
                                                         return_samples=True,
                                                         return_components=True)
# Remove systemic velocity
rv_model -= components['mu']
median_model -= components['mu']

insts = ['CORALIE','HARPS']
colors = ['royalblue','firebrick',]
symbols = ['.','^']
msizes = [14,10]

fig, (ax1, ax2) = plt.subplots(2, 1,sharex=True,
     gridspec_kw={'height_ratios': [2.5, 1]},figsize=(10,8))

for i in range(len(insts)):

    #RVS
    mu = np.median(results.posteriors['posterior_samples']['mu_'+insts[i]])
    jitter = np.median(results.posteriors['posterior_samples']['sigma_w_'+insts[i]])
    markers, caps, bars = ax1.errorbar(dataset.times_rv[insts[i]], (dataset.data_rv[insts[i]]-mu)*1000.0,
                                       yerr=dataset.errors_rv[insts[i]]*1000.0,color=colors[i],
                                       fmt=symbols[i], ms=msizes[i], elinewidth=4, label=insts[i])
    [bar.set_alpha(0.4) for bar in bars]
    
    # Plotting RV models
    ax1.plot(model_rv_times, median_model*1000.0, color='black',linewidth=2)
    
    #residuals
    residuals = (dataset.data_rv[insts[i]] - results.rv.evaluate(insts[i]))*1000.0
    markers, caps, bars = ax2.errorbar(dataset.times_rv[insts[i]], residuals,color=colors[i],
                                       yerr=dataset.errors_rv[insts[i]]*1000.0,
                                       fmt=symbols[i], ms=msizes[i], elinewidth=4, label=insts[i])
    [bar.set_alpha(0.4) for bar in bars]

ax2.axhline(0, color='black',alpha=1,linewidth=2)
ax2.set_ylim([-125,125])

ax1.set_ylabel('RV (m/s)',fontsize=28)
ax2.set_ylabel('residuals',fontsize=28)
ax2.set_xlabel('BJD - 2400000',fontsize=28)
ax1.legend(fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=24,width=3,length=10)
ax2.tick_params(axis='both', which='major', labelsize=24,width=3,length=10)

fig.tight_layout()
plt.savefig(output_folder+'/RVfull_'+star+'.pdf')

########################################################
############# PLOT RV periodogram  #####################
########################################################

from astropy.stats import LombScargle
from scipy import signal

#combine RVs with mu subtracted
rv_coralie = dataset.data_rv['CORALIE'] - np.median(results.posteriors['posterior_samples']['mu_CORALIE'])
rv_harps = dataset.data_rv['HARPS'] - np.median(results.posteriors['posterior_samples']['mu_HARPS'])
rv_offseted = np.concatenate([rv_coralie,rv_harps])
rverr_offseted = np.concatenate([dataset.errors_rv['CORALIE'],dataset.errors_rv['HARPS']])
rvtime_offseted = np.concatenate([dataset.times_rv['CORALIE'],dataset.times_rv['HARPS']])
#plt.errorbar(rvtime_offseted,rv_offseted,yerr=rverr_offseted,fmt='.')

max_period = 99999.9
min_period = 0.09
min_frequency = 1./max_period
max_frequency = 1./min_period

ls = LombScargle(rvtime_offseted,rv_offseted,rverr_offseted)
freq, power = ls.autopower(minimum_frequency=min_frequency, maximum_frequency=max_frequency)
period = 1./freq
#print(power.max())


#FAP
ls.false_alarm_probability(power.max())  
probabilities = [0.1, 0.01, 0.001]
faps = ls.false_alarm_level(probabilities)

plt.figure(figsize=(10,6))
plt.xscale('log')
plt.xlim(0.9,1100)
plt.ylim(0,1)

plt.axhline(y=faps[2], color='firebrick',linewidth=3, linestyle='--',label='0.1% FAP')
plt.legend(fontsize=20)
plt.title('RV periodogram',fontsize=32)
plt.plot(period, power,color='cornflowerblue',linewidth=4)
plt.xlabel('Period (days)',fontsize=28)
plt.ylabel('Normalized power',fontsize=28)
plt.xticks([1,10,100,1000],['1','10','100','1000'],fontsize=24)
plt.yticks(fontsize=24)
plt.axvline(np.median(P),linewidth=4,color='black',alpha=0.5,label=str(np.median(P).round(2)))
plt.text(np.median(P)-70,0.9,'P = '+str(np.median(P).round(2)),fontsize=24)
plt.tick_params(axis='both',width=3,length=10)
plt.tight_layout()
plt.savefig(output_folder+'/RV_periodogram_'+star+'.pdf')

#fap_1 = plt.axhline(y=faps[0], color='r', linestyle='--',label='10% FAP')
#fap_2 = plt.axhline(y=faps[1], color='g', linestyle='--',label='1% FAP')
#fap_3 =  plt.axhline(y=faps[2], color='b', linestyle='--',label='0.1% FAP')
#plt.legend(handles=[fap_1,fap_2,fap_3])

####### RESIDUALS PERIODOGRAM ######
residuals_coralie = (dataset.data_rv['CORALIE'] - results.rv.evaluate('CORALIE'))
residuals_harps = (dataset.data_rv['HARPS'] - results.rv.evaluate('HARPS'))
resrv_offseted = np.concatenate([residuals_coralie,residuals_harps])

ls = LombScargle(rvtime_offseted,resrv_offseted,rverr_offseted)
freq, power = ls.autopower(minimum_frequency=min_frequency, maximum_frequency=max_frequency)
period = 1./freq

#FAP
ls.false_alarm_probability(power.max())  
probabilities = [0.1, 0.01, 0.001]
faps = ls.false_alarm_level(probabilities)

plt.figure(figsize=(10,6))
plt.title('RV residuals periodogram',fontsize=32)
plt.xscale('log')
plt.xlim(0.9,1100)
plt.ylim(0,1)
plt.axhline(y=faps[2], color='firebrick',linewidth=3, linestyle='--',label='0.1% FAP')
plt.legend(fontsize=20)
plt.plot(period, power,color='cornflowerblue',linewidth=4)
plt.xlabel('Period (days)',fontsize=28)
plt.ylabel('Normalized power',fontsize=28)
plt.xticks([1,10,100,1000],['1','10','100','1000'],fontsize=24)
plt.yticks(fontsize=24)
plt.tick_params(axis='both',width=3,length=10)
plt.tight_layout()
plt.savefig(output_folder+'/RVresid_periodogram_'+star+'.pdf')

############# PHASED RVS with residuals #####################

insts = ['CORALIE','HARPS']
colors = ['royalblue','firebrick',]
symbols = ['.','^']
msizes = [14,10]


fig, (ax1, ax2) = plt.subplots(2, 1,sharex=True,
     gridspec_kw={'height_ratios': [2.5, 1]},figsize=(10,8))
    
for i in range(len(insts)):
    mu = np.median(results.posteriors['posterior_samples']['mu_'+insts[i]])
    jitter = np.median(results.posteriors['posterior_samples']['sigma_w_'+insts[i]])
    phases = juliet.get_phases(dataset.times_rv[insts[i]], P, t0)
    markers, caps, bars = ax1.errorbar(phases, (dataset.data_rv[insts[i]]-mu)*1000.0,
                                       yerr=(dataset.errors_rv[insts[i]])*1000.0,color=colors[i],
                                       fmt=symbols[i], ms=msizes[i], elinewidth=4, label=insts[i])
    [bar.set_alpha(0.5) for bar in bars]

    #residuals
    residuals = (dataset.data_rv[insts[i]] - results.rv.evaluate(insts[i]))*1000.0
    markers, caps, bars = ax2.errorbar(phases, residuals,color=colors[i],
                                       yerr=dataset.errors_rv[insts[i]]*1000.0,
                                       fmt=symbols[i],ms=msizes[i],elinewidth=4, label=insts[i])
    [bar.set_alpha(0.5) for bar in bars]

ax2.axhline(0, color='black',alpha=1,linewidth=2)

# Plotting RV models
phases = juliet.get_phases(model_rv_times, P, t0)
idx = np.argsort(phases)
ax1.plot(phases[idx], median_model[idx]*1000.0, color='black',linewidth=2)

ax1.text(0.2,120,'P = %.1f days'%per_value,fontsize=22)#,transform = plt.transAxes)
ax1.text(0.2,90,'K = %.1f m/s'%K_value,fontsize=22)#,transform = plt.transAxes)
ax1.text(0.2,60,'ecc = %.2f'%ecc_value,fontsize=22)#,transform = plt.transAxes)

ax1.set_ylabel('RV (m/s)',fontsize=28)
ax2.set_ylabel('residuals',fontsize=28)
ax2.set_xlabel('phase',fontsize=28)
ax1.legend(fontsize=20)
ax1.tick_params(axis='y', which='major', labelsize=24,width=3,length=10)
ax2.tick_params(axis='both', which='major', labelsize=24,width=3,length=10)

plt.tight_layout()
plt.savefig(output_folder+'/RVphase_withResiduals_'+star+'.pdf')

############# PHASED RVS #####################

insts = ['CORALIE','HARPS']
colors = ['royalblue','firebrick',]
symbols = ['.','^']
msizes = [14,10]


fig, (ax1) = plt.subplots(figsize=(10,6))
    
for i in range(len(insts)):
    mu = np.median(results.posteriors['posterior_samples']['mu_'+insts[i]])
    jitter = np.median(results.posteriors['posterior_samples']['sigma_w_'+insts[i]])
    phases = juliet.get_phases(dataset.times_rv[insts[i]], P, t0)
    markers, caps, bars = ax1.errorbar(phases, (dataset.data_rv[insts[i]]-mu)*1000.0,
                                       yerr=(dataset.errors_rv[insts[i]])*1000.0,color=colors[i],
                                       fmt=symbols[i], ms=msizes[i], elinewidth=4, label=insts[i])
    [bar.set_alpha(0.5) for bar in bars]


# Plotting RV models
phases = juliet.get_phases(model_rv_times, P, t0)
idx = np.argsort(phases)
ax1.plot(phases[idx], median_model[idx]*1000.0, color='black',linewidth=2)

ax1.text(-0.47,-130,'P = %.2f days'%per_value,fontsize=24)#,transform = plt.transAxes)
ax1.text(-0.47,-160,'K = %.1f m/s'%K_value,fontsize=24)#,transform = plt.transAxes)
ax1.text(-0.47,-190,'ecc = %.2f'%ecc_value,fontsize=24)#,transform = plt.transAxes)

ax1.set_ylabel('RV (m/s)',fontsize=28)
ax1.set_xlabel('phase',fontsize=28)
ax1.legend(fontsize=24)
ax1.tick_params(axis='y', which='major', labelsize=24,width=3,length=10)
ax1.tick_params(axis='both', which='major', labelsize=24,width=3,length=10)

plt.tight_layout()
plt.savefig(output_folder+'/RVphase_'+star+'.pdf')


########## COMBINE PERIODOGRAM########

fig, (ax1, ax2) = plt.subplots(2, 1,sharex=True,
     gridspec_kw={'height_ratios': [1, 1]},figsize=(11,8))

#combine RVs with mu subtracted
rv_coralie = dataset.data_rv['CORALIE'] - np.median(results.posteriors['posterior_samples']['mu_CORALIE'])
rv_harps = dataset.data_rv['HARPS'] - np.median(results.posteriors['posterior_samples']['mu_HARPS'])
rv_offseted = np.concatenate([rv_coralie,rv_harps])
rverr_offseted = np.concatenate([dataset.errors_rv['CORALIE'],dataset.errors_rv['HARPS']])
rvtime_offseted = np.concatenate([dataset.times_rv['CORALIE'],dataset.times_rv['HARPS']])
#plt.errorbar(rvtime_offseted,rv_offseted,yerr=rverr_offseted,fmt='.')

max_period = 99999.9
min_period = 0.09
min_frequency = 1./max_period
max_frequency = 1./min_period

ls = LombScargle(rvtime_offseted,rv_offseted,rverr_offseted)
freq, power = ls.autopower(minimum_frequency=min_frequency, maximum_frequency=max_frequency)
period = 1./freq


#FAP
ls.false_alarm_probability(power.max())  
probabilities = [0.1, 0.01, 0.001]
faps = ls.false_alarm_level(probabilities)

ax1.set_xscale('log')
ax1.set_xlim(0.9,1100)
ax1.set_ylim(0,1)

ax1.axhline(y=faps[2], color='firebrick',linewidth=3, linestyle='--',label='0.1% FAP')
ax1.legend(fontsize=22)
ax1.plot(period, power,color='cornflowerblue',linewidth=4)
ax1.axvline(np.median(P),linewidth=4,color='black',alpha=0.5,label='P = '+str(np.median(P).round(2)))


####### RESIDUALS PERIODOGRAM ######
residuals_coralie = (dataset.data_rv['CORALIE'] - results.rv.evaluate('CORALIE'))
residuals_harps = (dataset.data_rv['HARPS'] - results.rv.evaluate('HARPS'))
resrv_offseted = np.concatenate([residuals_coralie,residuals_harps])

ls = LombScargle(rvtime_offseted,resrv_offseted,rverr_offseted)
freq, power = ls.autopower(minimum_frequency=min_frequency, maximum_frequency=max_frequency)
period = 1./freq

#FAP
ls.false_alarm_probability(power.max())  
probabilities = [0.1, 0.01, 0.001]
faps = ls.false_alarm_level(probabilities)

#ax2.xscale('log')
ax2.set_ylim(0,1)
ax2.axhline(y=faps[2], color='firebrick',linewidth=3, linestyle='--',label='0.1% FAP')
ax2.legend(fontsize=22)
ax2.plot(period, power,color='cornflowerblue',linewidth=4)
ax2.set_xticklabels(['0.01','0.1','1','10','100','1000','10000','100000'])
fig.text(-0.0, 0.5, 'Normalized power', va='center', rotation='vertical',fontsize=30)

#ax2.set_ylabel('Normalized power',fontsize=28)
ax2.set_xlabel('Period (days)',fontsize=28)
ax1.legend(fontsize=22)
ax1.tick_params(axis='both', which='major', labelsize=24,width=3,length=10)
ax2.tick_params(axis='both', which='major', labelsize=24,width=3,length=10)

ax1.text(1.5,0.8,'RV periodogram',fontsize=28)
ax2.text(1.5,0.8,'RV residuals periodogram',fontsize=28)
#plt.tight_layout()
plt.savefig(output_folder+'/RVperiodograms_'+star+'.pdf')

########################################################
################       LIGHTCURVES      ################ 
########################################################

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
		tot = np.sum(x / y for x, y in zip(data[lims], error[lims])) / np.sum(1.0/error[lims]) 
		tot_error = sum(error[lims]**2)
		
		if n > 0:
			bmids.append(mid)
			bvals.append(tot)
			#berror.append(np.sqrt(tot_error)/n)
            ##### NOLAN CHANGE berror to standard deviation of points
			#berror.append(np.nanstd(data[lims]))
            #NG change again to the standard error
			berror.append(np.nanstd(data[lims])/np.sqrt(n))
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

########################################################
################ PLOT PHASED LIGHTCURVE ################ 
########################################################
# Extract new period and time-of-transit center:
P,t0 =  np.median(results.posteriors['posterior_samples']['P_p1']),\
        np.median(results.posteriors['posterior_samples']['t0_p1'])

# Generate arrays to super-sample the models:
model_phases = np.linspace(-0.04,0.04,1000)
model_times = model_phases*P + t0

# Plot figure:
#fig, ax1  = plt.subplots(figsize=(10, 18))

fig, (ax1, ax2) = plt.subplots(2, 1,sharex=True,
     gridspec_kw={'height_ratios': [3, 1]},figsize=(10,18))

##### TESS13 #####
instrument = 'TESS13'
offset = 0
model_lc_full = results.lc.evaluate(instrument)
model_lc = results.lc.model[instrument]['deterministic']
model_lc_gp = results.lc.model[instrument]['GP']
phases = juliet.utils.get_phases(dataset.times_lc[instrument], P, t0)

residuals_tess13 = dataset.data_lc[instrument] - model_lc_full

ax1.scatter(phases*P*24.0, dataset.data_lc[instrument]-offset-model_lc_gp,color='dimgrey',s=20,edgecolor='black',alpha=1,label=instrument)
ax1.errorbar(phases*P*24.0, dataset.data_lc[instrument]-offset-model_lc_gp, \
            yerr = dataset.errors_lc[instrument], alpha = 0.5,markersize='3',color='black',fmt='o',linewidth=1)
ax1.plot(phases*P*24.0,model_lc-offset, color='firebrick',zorder=10,linewidth=2,label='GP + Transit Model')
ax1.text(-9.5,1.004,r'TESS13 30-min 14 July 2019',fontsize=26)

##### TESS27 #####
instrument = 'TESS27'
offset = 0.03
model_lc_full = results.lc.evaluate(instrument)
model_lc = results.lc.model[instrument]['deterministic']
model_lc_gp = results.lc.model[instrument]['GP']
phases = juliet.utils.get_phases(dataset.times_lc[instrument], P, t0)

ax1.scatter(phases*P*24.0, dataset.data_lc[instrument]-offset-model_lc_gp,color='dimgrey',s=20,edgecolor='black',alpha=1,label=instrument)
ax1.errorbar(phases*P*24.0, dataset.data_lc[instrument]-offset-model_lc_gp, \
            yerr = dataset.errors_lc[instrument], alpha = 0.5,markersize='3',color='black',fmt='o',linewidth=1)
ax1.plot(phases*P*24.0,model_lc-offset, color='firebrick',zorder=10,linewidth=2,label='GP + Transit Model')
ax1.text(-9.5,0.982,'TESS27 2-min 24 July 2020',fontsize=26)

##### SAAO #####
binwidth = 2.0/60.0/24.0 #in days
saao_bin_time, saao_bin_flux, saao_bin_fluxerr = fancy_bin(dataset.times_lc['SAAO'],dataset.data_lc['SAAO'],dataset.errors_lc['SAAO'],binwidth)
instrument = 'SAAO'
offset = 0.06
model_lc = results.lc.evaluate(instrument, t = model_times)
#phases = juliet.utils.get_phases(dataset.times_lc[instrument], P, t0)
phases = juliet.utils.get_phases(saao_bin_time, P, t0)
ax1.scatter(phases*P*24.0, saao_bin_flux-offset,color='dimgrey',s=20,edgecolor='black',alpha=1,label=instrument)
ax1.errorbar(phases*P*24.0, saao_bin_flux-offset, \
            yerr = saao_bin_fluxerr, alpha = 0.5,markersize='3',color='black',fmt='o',linewidth=1)
ax1.plot(model_phases*P*24.0,model_lc-offset, color='firebrick',zorder=10,linewidth=2,label=' Model')
ax1.text(-9.5,0.9465,r'SAAO V 2-min bin 3 August 2021',fontsize=26)

##### ECAM #####
instrument = 'ECAM'
offset = 0.09
model_lc = results.lc.evaluate(instrument, t = model_times)
phases = juliet.utils.get_phases(dataset.times_lc[instrument], P, t0)

binwidth = 30.0/60.0/24.0 #in days
ecam_bin_time,ecam_bin_flux, ecam_bin_fluxerr = fancy_bin(dataset.times_lc['ECAM'],dataset.data_lc['ECAM'],dataset.errors_lc['ECAM'],binwidth)
phases_ecam_bin = juliet.utils.get_phases(ecam_bin_time, P, t0)

ax1.scatter(phases*P*24.0, dataset.data_lc[instrument]-offset,color='dimgrey',s=20,edgecolor='black',alpha=1,label=instrument)
ax1.errorbar(phases*P*24.0, dataset.data_lc[instrument]-offset, \
            yerr = dataset.errors_lc[instrument], alpha = 0.5,markersize='3',color='black',fmt='o',linewidth=1)
ax1.plot(model_phases*P*24.0,model_lc-offset, color='firebrick',zorder=10,linewidth=2,label='Transit Model')
ax1.text(-9.5,0.921,r'ECAM V 2.5-min 3 August 2021',fontsize=26)

###################
#### residuals ####
###################
instrument = 'TESS13'
offset=0.0
model_lc = results.lc.evaluate(instrument)
phases = juliet.utils.get_phases(dataset.times_lc[instrument], P, t0)
ax2.scatter(phases*P*24.0, dataset.data_lc[instrument]-model_lc+offset,
            color='dimgrey',s=20,edgecolor='black',alpha=1,label=instrument)
ax2.errorbar(phases*P*24.0, dataset.data_lc[instrument]-model_lc+offset,
            yerr = dataset.errors_lc[instrument], alpha = 0.5,markersize='3',
            color='black',fmt='o',linewidth=1)
ax2.axhline(offset,linestyle='--',color='firebrick',zorder=10,linewidth=2,label='')

resid = dataset.data_lc[instrument]-model_lc
rms = np.sqrt(np.mean(resid**2))
ax2.text(-9.5,offset+0.0065,r'TESS13 30-min 14 July 2019 (RMS = %.5f)'%(rms),fontsize=18,color='black',zorder=100)

#### residuals ####
instrument = 'TESS27'
offset=-0.03
model_lc = results.lc.evaluate(instrument)
phases = juliet.utils.get_phases(dataset.times_lc[instrument], P, t0)
ax2.scatter(phases*P*24.0, dataset.data_lc[instrument]-model_lc+offset,
            color='dimgrey',s=20,edgecolor='black',alpha=1,label=instrument)
ax2.errorbar(phases*P*24.0, dataset.data_lc[instrument]-model_lc_full+offset,
            yerr = dataset.errors_lc[instrument], alpha = 0.5,markersize='3',
            color='black',fmt='o',linewidth=1)
ax2.axhline(offset,linestyle='--', color='firebrick',zorder=10,linewidth=2,label='')

#resid = dataset.data_lc[instrument]-model_lc_full
#rms = np.sqrt(np.mean(resid**2))
binwidth = 30.0/60.0/24.0 #in days
instrument = 'TESS27'
bin_time, bin_flux, bin_fluxerr = fancy_bin(dataset.times_lc[instrument],dataset.data_lc[instrument],dataset.errors_lc[instrument],binwidth)
resid = bin_flux -  results.lc.evaluate(instrument, t = bin_time, GPregressors = bin_time)
rms = np.sqrt(np.mean(resid**2))
ax2.text(-9.5,offset+0.013,'TESS27 2-min 24 July 2020 (30-min bin RMS = %.5f)'%(rms),fontsize=18,color='black',zorder=100)

#### residuals ####
instrument = 'SAAO'
offset=-0.06
model_lc = results.lc.evaluate(instrument, t = saao_bin_time)
phases = juliet.utils.get_phases(saao_bin_time, P, t0)
ax2.scatter(phases*P*24.0, saao_bin_flux-model_lc+offset,color='dimgrey',s=20,edgecolor='black',alpha=1,label=instrument)
ax2.errorbar(phases*P*24.0, saao_bin_flux-model_lc+offset,yerr=saao_bin_fluxerr,
             alpha=0.5, markersize='3', color='black',fmt='o',linewidth=1)
ax2.axhline(offset,linestyle='--', color='firebrick',zorder=10,linewidth=2,label=' Model')


binwidth = 30.0/60.0/24.0 #in days
instrument = 'SAAO'
bin_time, bin_flux, bin_fluxerr = fancy_bin(dataset.times_lc[instrument],dataset.data_lc[instrument],dataset.errors_lc[instrument],binwidth)
resid = bin_flux -  results.lc.evaluate(instrument, t = bin_time)
rms = np.sqrt(np.mean(resid**2))
ax2.text(-9.5,offset+0.0075,r'SAAO V 2-min bin 3 August 2021 (30-min bin RMS = %.5f)'%(rms),fontsize=18,color='black',zorder=100)

#### residuals ####
instrument = 'ECAM'
offset=-0.09
model_lc = results.lc.evaluate(instrument)
phases = juliet.utils.get_phases(dataset.times_lc[instrument], P, t0)

ax2.scatter(phases*P*24.0, dataset.data_lc[instrument]-model_lc+offset,color='dimgrey',
            s=20,edgecolor='black',alpha=1,label=instrument)
ax2.errorbar(phases*P*24.0, dataset.data_lc[instrument]-model_lc+offset,yerr = dataset.errors_lc[instrument],
             alpha = 0.5,markersize='3',color='black',fmt='o',linewidth=1)
ax2.axhline(offset,linestyle='--', color='firebrick',zorder=10,linewidth=2,label='')

#resid = dataset.data_lc[instrument]-model_lc
#rms = np.sqrt(np.mean(resid**2))
binwidth = 30.0/60.0/24.0 #in days
instrument = 'ECAM'
bin_time, bin_flux, bin_fluxerr = fancy_bin(dataset.times_lc[instrument],dataset.data_lc[instrument],dataset.errors_lc[instrument],binwidth)
resid = bin_flux - results.lc.evaluate(instrument, t = bin_time)
rms = np.sqrt(np.mean(resid**2))
ax2.text(-9.5,offset+0.0075,r'ECAM V 2.5-min 3 August 2021 (30-min bin RMS = %.5f)'%(rms),fontsize=18,color='black',zorder=100)

#####
ax1.tick_params(axis='both',width=3,length=10)
plt.xticks(fontsize=27)
ax1.set_yticks([0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.0])
ax1.set_yticklabels(['0.90','0.91','0.92','0.93','0.94','0.95','0.96','0.97','0.98','0.99','1.00'],fontsize=27)
ax2.set_yticks([-0.10,-0.08,-0.06,-0.04,-0.02,0.00,0.02])
ax2.set_yticklabels(['-0.10','-0.08','-0.06','-0.04','-0.02','0.00','0.02'],fontsize=27)
plt.xlim(-10,10)
plt.xlabel(r'Time - T$_{0}$ (hours)',fontsize=30)
ax1.set_ylabel('Relative flux + offset',fontsize=30)
ax2.set_ylabel('Residuals + offset',fontsize=30)
#plt.legend(loc='lower left',fontsize=27,framealpha=1,facecolor='white')#,fancybox=True,shadow=True)
plt.tick_params(axis='both',width=3,length=10)
plt.tight_layout()
plt.savefig(output_folder+'/transit_phase_'+star+'.pdf')


########################################################
########################################################
########################################################
### Full Light Curves ###

# ### SAAO TRANSIT ###
# t0 = np.median(results.posteriors['posterior_samples']['t0_p1'])
# instrument='SAAO'
# # Plot. First extract model:
# transit_model, transit_up68, transit_low68, components  = results.lc.evaluate(instrument, return_err=True, \
#                                                                               return_components = True, \
#                                                                               all_samples = True)

# timeoffset = 59430

# plt.figure(figsize=(12,8))

# #plt.title(instrument,fontsize=24)
# plt.errorbar(dataset.times_lc[instrument]-timeoffset, dataset.data_lc[instrument], \
#              yerr = dataset.errors_lc[instrument], fmt = '.' , linewidth=0.5,alpha = 0.5,color='black',label=instrument)

# # Out-of-transit flux:
# oot_flux = np.median(1./(1. + results.posteriors['posterior_samples']['mflux_'+instrument]))

# # Plot non-transit model::
# plt.plot(dataset.times_lc[instrument]-timeoffset, oot_flux + components['lm'], color='cornflowerblue', lw = 3, label='detrending model')#label = 'Linear model + oot flux')
# plt.plot(dataset.times_lc[instrument]-timeoffset, transit_model, color='red', label = 'Full model (transit+detrending)')
# plt.fill_between(dataset.times_lc[instrument]-timeoffset,transit_up68,transit_low68,\
#                  color='red',alpha=0.5,zorder=5)

# plt.xlabel('BJD - 2459430',fontsize=28)
# plt.ylabel('Relative flux',fontsize=28)
# plt.xticks(fontsize=24)
# plt.yticks(fontsize=24)
# plt.tick_params(axis='both',width=3,length=10)
# plt.legend(fontsize=24)
# plt.savefig(output_folder+'/transit_detrend_'+instrument+'_'+star+'.pdf')

# instrument='ECAM'

# transit_model, transit_up68, transit_low68, components  = results.lc.evaluate(instrument, return_err=True, \
#                                                                               return_components = True, \
#                                                                               all_samples = True)

# timeoffset = 59430

# plt.figure(figsize=(12,8))

# #plt.title(instrument,fontsize=24)
# plt.errorbar(dataset.times_lc[instrument]-timeoffset, dataset.data_lc[instrument], \
#              yerr = dataset.errors_lc[instrument], fmt = '.' , linewidth=0.5,alpha = 0.5,color='black',label=instrument)

# # Out-of-transit flux:
# oot_flux = np.median(1./(1. + results.posteriors['posterior_samples']['mflux_'+instrument]))

# # Plot non-transit model::
# plt.plot(dataset.times_lc[instrument]-timeoffset, oot_flux + components['lm'], color='cornflowerblue', lw = 3, label='detrending model')#label = 'Linear model + oot flux')
# plt.plot(dataset.times_lc[instrument]-timeoffset, transit_model, color='red', label = 'Full model (transit+detrending)')
# plt.fill_between(dataset.times_lc[instrument]-timeoffset,transit_up68,transit_low68,\
#                  color='red',alpha=0.5,zorder=5)

# plt.xlabel('BJD - 2459430',fontsize=28)
# plt.ylabel('Relative flux',fontsize=28)
# plt.xticks(fontsize=24)
# plt.yticks(fontsize=24)
# plt.tick_params(axis='both',width=3,length=10)
# plt.legend(fontsize=24)

# plt.savefig(output_folder+'/transit_detrend_'+instrument+'_'+star+'.pdf')


########################################################



######## TESS 27 GP full PLOT ######## 


instrument = 'TESS27'
#transit_model, transit_up68, transit_low68  = results.lc.evaluate(instrument, GPregressors=dataset.times_lc[instrument])

transit_model = results.lc.evaluate(instrument,t=dataset.times_lc[instrument],GPregressors=dataset.times_lc[instrument])

fig = plt.figure(figsize=(16,8))
gs = gridspec.GridSpec(1, 2, width_ratios=[2,1])

ax1 = plt.subplot(gs[0])
ax1.scatter(dataset.times_lc[instrument], dataset.data_lc[instrument],color='dimgrey',s=10,edgecolor='black',alpha=0.9,label=instrument)
ax1.errorbar(dataset.times_lc[instrument], dataset.data_lc[instrument], \
            yerr = dataset.errors_lc[instrument], alpha = 0.1,markersize='3',color='black',fmt='o',linewidth=1)
ax1.plot(dataset.times_lc[instrument],transit_model, color='firebrick',zorder=100,linewidth=3,label='GP + Transit Model')


T0tess = 59054.96842


ax2 = plt.subplot(gs[1])
ax2.scatter(dataset.times_lc[instrument]-T0tess, dataset.data_lc[instrument],color='dimgrey',s=10,edgecolor='black',alpha=0.9,label=instrument)
ax2.errorbar(dataset.times_lc[instrument]-T0tess, dataset.data_lc[instrument], \
            yerr = dataset.errors_lc[instrument], alpha = 0.1,markersize='3',color='black',fmt='o',linewidth=1)
ax2.plot(dataset.times_lc[instrument]-T0tess,transit_model, color='firebrick',zorder=100,linewidth=3,label='GP + Transit Model')


ax2.yaxis.set_major_formatter(plt.NullFormatter())
ax2.set_xlim([-1.5,1.5])

ax1.tick_params(which='both', width=1.5,labelsize=24)
ax1.tick_params(which='major', length=12)
ax1.tick_params(which='minor', length=8)#, color='r')
ax1.tick_params(which='both', direction='out')

ax2.tick_params(which='both', width=1.5,labelsize=24)
ax2.tick_params(which='major', length=12)
ax2.tick_params(which='minor', length=8)#, color='r')
ax2.tick_params(which='both', direction='out')

ax1.set_xlabel('BJD - 2400000',fontsize=28)
ax2.set_xlabel(r'Time - T$_{0}$ (days)',fontsize=28)
ax1.set_ylabel('Relative flux',fontsize=28)
ax1.legend(loc='lower left',fontsize=24)
plt.tight_layout()

plt.savefig(output_folder+'/transit_full_'+star+'_TESS27.pdf')


######## TESS 13 GP Full PLOT ######## 
instrument = 'TESS13'
#transit_model, transit_up68, transit_low68  = results.lc.evaluate(instrument, t = dataset.times_lc[instrument],return_err=True)
#transit_model = results.lc.evaluate(instrument)
transit_model = results.lc.evaluate(instrument,t=dataset.times_lc[instrument],GPregressors=dataset.times_lc[instrument])


fig = plt.figure(figsize=(16,8))
gs = gridspec.GridSpec(1, 2, width_ratios=[2,1])

ax1 = plt.subplot(gs[0])
ax1.scatter(dataset.times_lc[instrument], dataset.data_lc[instrument],color='dimgrey',s=10,edgecolor='black',alpha=0.9,label=instrument)
ax1.errorbar(dataset.times_lc[instrument], dataset.data_lc[instrument], \
            yerr = dataset.errors_lc[instrument], alpha = 0.1,markersize='3',color='black',fmt='o',linewidth=1)
ax1.plot(dataset.times_lc[instrument],transit_model, color='firebrick',zorder=100,linewidth=3,label='GP + Transit Model')


T0tess = 58679.35107


ax2 = plt.subplot(gs[1])
ax2.scatter(dataset.times_lc[instrument]-T0tess, dataset.data_lc[instrument],color='dimgrey',s=10,edgecolor='black',alpha=0.9,label=instrument)
ax2.errorbar(dataset.times_lc[instrument]-T0tess, dataset.data_lc[instrument], \
            yerr = dataset.errors_lc[instrument], alpha = 0.1,markersize='3',color='black',fmt='o',linewidth=1)
ax2.plot(dataset.times_lc[instrument]-T0tess,transit_model, color='firebrick',zorder=100,linewidth=3,label='GP + Transit Model')


ax2.yaxis.set_major_formatter(plt.NullFormatter())
ax2.set_xlim([-1.5,1.5])

ax1.tick_params(which='both', width=1.5,labelsize=24)
ax1.tick_params(which='major', length=12)
ax1.tick_params(which='minor', length=8)#, color='r')
ax1.tick_params(which='both', direction='out')

ax2.tick_params(which='both', width=1.5,labelsize=24)
ax2.tick_params(which='major', length=12)
ax2.tick_params(which='minor', length=8)#, color='r')
ax2.tick_params(which='both', direction='out')

ax1.set_xlabel('BJD - 2400000',fontsize=28)
ax2.set_xlabel(r'Time - T$_{0}$ (days)',fontsize=28)
ax1.set_ylabel('Relative flux',fontsize=28)
ax1.legend(loc='lower left',fontsize=24)
plt.tight_layout()

plt.savefig(output_folder+'/transit_full_'+star+'_TESS13.pdf')