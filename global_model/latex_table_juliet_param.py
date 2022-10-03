#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create latex table for Juliet global model parameters
Created on Fri Sep  2 15:39:25 2022
@author: nolangrieves
"""

import pandas as pd
import numpy as np
from juliet_definitions import *
import pykima as pk

outputdir = 'juliet_output_TOI-5542_1sep2022'

####### Stellar Parameters #######
Mstar = 0.890
Mstar_err = 0.057
Rstar = 1.058 #* R_sun
Rstar_err = 0.036 #* R_sun #-0.019 +0.030

Lstar = 1.057
Lstar_err = 0.0593

Teff = 5700
Teff_err = 80


#################################################################
#### Load Posteriors Output file
#################################################################

post = pd.read_csv(outputdir+'/posteriors.dat',skiprows=[0],header=None,sep=r"\s+",index_col=0)
post = post.T

#################################################################
#### Fit/Load Actual Posteriors
#################################################################
dataset = juliet.load(priors=outputdir+'/priors.dat',
                      lcfilename=outputdir+'/lc.dat',
                      GPlceparamfile=outputdir+'/GP_lc_regressors.dat',
                      rvfilename=outputdir+'/rvs.dat',
                      out_folder=outputdir)

results = dataset.fit(n_live_points=400,sampler='dynesty',nthreads=4)

#################################################################
#### Derive Results from Actual Posteriors ###
#################################################################

### Planet Radius ###
r1 = results.posteriors['posterior_samples']['r1_p1']  
r2 = results.posteriors['posterior_samples']['r2_p1']  
b,p = juliet.utils.reverse_bp(r1, r2, 0., 1.) # impact parameter b = (a / Rs) *cos(i) # planet to star ratio k = (Rp / Rs)
#b = results.posteriors['posterior_samples']['b_p1']  # impact parameter b = (a / Rs) *cos(i)
#p = results.posteriors['posterior_samples']['p_p1']  # planet to star ratio k = (Rp / Rs)

Rstar_array = np.random.normal(Rstar,Rstar_err,len(p)) 
radius_array = p*Rstar_array * 9.73116 #Rsun to Rjupiter 109.076 #Rsun to Rearth  # planet to star ratio k = (Rp / Rs)

Rp_value = get_quantiles(radius_array)[0]
Rp_value_up = (get_quantiles(radius_array)[1]) - (get_quantiles(radius_array)[0])
Rp_value_low = abs((get_quantiles(radius_array)[2]) - (get_quantiles(radius_array)[0]))

### Planet Mass ###
# get inclination from eccentricity and omega from Juliet
ecc = results.posteriors['posterior_samples']['ecc_p1']
omega = results.posteriors['posterior_samples']['omega_p1']

# get scaled semi major axis from stellar_density
# from juliet utils.py, function writepp
G = 6.67408e-11
try:
    a = ((results.posteriors['posterior_samples']['rho']
          *G*((results.posteriors['posterior_samples']['P_'+'p1']*24.*3600.)**2))/(3.*np.pi))**(1./3.)
except KeyError:
    a = results.posteriors['posterior_samples']['a_p1']

ecc_factor = (1. + ecc*np.sin(omega))/(1. - ecc**2)
inc_inv_factor = (b/a)*ecc_factor
inc = np.arccos(inc_inv_factor)*180./np.pi

Mpsini = pk.utils.get_planet_mass(results.posteriors['posterior_samples']['P_p1'],
                                  results.posteriors['posterior_samples']['K_p1']*1000.0,
                                  e=results.posteriors['posterior_samples']['ecc_p1'],
                                  star_mass=np.random.normal(Mstar,Mstar_err,len(results.posteriors['posterior_samples']['P_p1'])), full_output=True)
Mp = Mpsini[2] / np.sin((inc*np.pi/180.))
Mp_value = get_quantiles(Mp)[0]
Mp_value_up = (get_quantiles(Mp)[1]) - (get_quantiles(Mp)[0])
Mp_value_low = abs((get_quantiles(Mp)[2]) - (get_quantiles(Mp)[0]))

#Jupiter radius = 7.1492×10^7 m
#Jupiter mass = 1.89813x10^27 kg
planet_density = (Mp*1.89813e30) / ((4.0/3.0)*np.pi*(radius_array * 7.1492e9)**3.0) #Jupiter radius = 7.1492×107 m
planet_density_value = get_quantiles(planet_density)[0]
planet_density_value_up = (get_quantiles(planet_density)[1]) - (get_quantiles(planet_density)[0])
planet_density_value_low = abs((get_quantiles(planet_density)[2]) - (get_quantiles(planet_density)[0]))

# a_p1 = a/R*

### a_AU = a*Rstar_array #get a in AU
#1 Rsun = 0.0046524726 AU
a_AU = a*Rstar_array*0.0046524726 #get a in AU
a_AU_value = get_quantiles(a_AU)[0]
a_AU_value_up = (get_quantiles(a_AU)[1]) - (get_quantiles(a_AU)[0])
a_AU_value_low = abs((get_quantiles(a_AU)[2]) - (get_quantiles(a_AU)[0]))

### Planet Insolation ###
#https://exoplanetarchive.ipac.caltech.edu/docs/poet_calculations.html
Lstar_array = np.random.normal(Lstar,Lstar_err,len(a))


S = Lstar / a_AU**2.0
S_value = get_quantiles(S)[0]
S_value_up = (get_quantiles(S)[1]) - (get_quantiles(S)[0])
S_value_low = abs((get_quantiles(S)[2]) - (get_quantiles(S)[0]))

### Equilibirium Temperature ###
#https://en.wikipedia.org/wiki/Planetary_equilibrium_temperature
bond_albedo = 0.343 #https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html
#1 Rsun = 0.0046524726 AU

Tstar_array =  np.random.normal(Teff,Teff_err,len(a))

Teq = Tstar_array*np.sqrt(Rstar_array*0.0046524726/(2.0*a_AU))*(1 - bond_albedo)**(1.0/4.0)
Teq_value = get_quantiles(Teq)[0]
#print(get_quantiles(Teq)[0])
#Teq_value_up = (get_quantiles(Teq)[1]) - (get_quantiles(Teq)[0])
#Teq_value_low = abs((get_quantiles(Teq)[2]) - (get_quantiles(Teq)[0]))

#calculate using 0, 0.343, 0.686
bond_albedo = 0
Teq = Tstar_array*np.sqrt(Rstar_array*0.0046524726/(2.0*a_AU))*(1 - bond_albedo)**(1.0/4.0)
Teq_value_up = get_quantiles(Teq)[0]-Teq_value
#print(get_quantiles(Teq)[0],get_quantiles(Teq)[0]-Teq_value)

bond_albedo = 0.686
Teq = Tstar_array*np.sqrt(Rstar_array*0.0046524726/(2.0*a_AU))*(1 - bond_albedo)**(1.0/4.0)
Teq_value_low = Teq_value-get_quantiles(Teq)[0]
#print(get_quantiles(Teq)[0],Teq_value-get_quantiles(Teq)[0])





fout=open(outputdir+'/TOI-5542_juliet_param.tex', 'w')
fout.write(r'\begin{table*}'+ "\n")
fout.write(r''+"\n")

fout.write(r'~~~~$q_{\text{1,TESS}}$\dotfill & Quadratic limb-darkening parametrization\dotfill & $\mathcal{N}(0.335,0.011)$\dotfill & $%.3f^{+%.3f}_{-%.3f}$ \\'
           %(post['q1_TESS13_TESS27'][1],post['q1_TESS13_TESS27'][2],post['q1_TESS13_TESS27'][3])+"\n")
fout.write(r'~~~~$q_{\text{2,TESS}}$\dotfill & Quadratic limb-darkening parametrization\dotfill & $\mathcal{N}(0.261,0.030)$\dotfill & $%.3f^{+%.3f}_{-%.3f}$ \\'
           %(post['q2_TESS13_TESS27'][1],post['q2_TESS13_TESS27'][2],post['q2_TESS13_TESS27'][3])+"\n")


fout.write(r'~~~~$m_{\text{flux,TESS13}}$\dotfill & Offset (relative flux)\dotfill & $\mathcal{N}(0,0.01)$\dotfill & $%.5f^{+%.5f}_{-%.5f}$ \\'
           %(post['mflux_TESS13'][1],post['mflux_TESS13'][2],post['mflux_TESS13'][3])+"\n")
#fout.write(r'~~~~$\sigma_{\text{TESS13}}$\dotfill  & Jitter (ppm)\dotfill & $\mathcal{J}(10^{-5},100)$\dotfill & $%.2f^{+%.2f}_{-%.2f}$ \\'
#           %(post['sigma_w_TESS13'][1],post['sigma_w_TESS13'][2],post['sigma_w_TESS13'][3])+"\n")

fout.write(r'~~~~$\sigma_{\text{GP,TESS13}}$\dotfill & GP amplitude (relative flux)\dotfill & $\mathcal{J}(10^{-6},1)$\dotfill &  $%.5f^{+%.5f}_{-%.5f}$ \\'
           %(post['GP_sigma_TESS13'][1],post['GP_sigma_TESS13'][2],post['GP_sigma_TESS13'][3])+"\n")
fout.write(r'~~~~$\rho_{\text{GP,TESS13}}$\dotfill & GP time-scale (days)\dotfill & $\mathcal{J}(10^{-6},10^{3})$\dotfill &  $%.2f^{+%.2f}_{-%.2f}$ \\'
           %(post['GP_rho_TESS13'][1],post['GP_rho_TESS13'][2],post['GP_rho_TESS13'][3])+"\n")

fout.write(r'~~~~$m_{\text{flux,TESS27}}$\dotfill & Offset (relative flux)\dotfill & $\mathcal{N}(0,0.01)$\dotfill &  $%.5f^{+%.5f}_{-%.5f}$ \\'
           %(post['mflux_TESS27'][1],post['mflux_TESS27'][2],post['mflux_TESS27'][3])+"\n")
#fout.write(r'~~~~$\sigma_{\text{TESS27}}$\dotfill  & Jitter (ppm)\dotfill & $\mathcal{J}(10^{-5},100)$\dotfill & $%.2f^{+%.2f}_{-%.2f}$ \\' 
#           %(post['sigma_w_TESS27'][1],post['sigma_w_TESS27'][2],post['sigma_w_TESS27'][3])+"\n")


fout.write(r'~~~~$\sigma_{\text{GP,TESS27}}$\dotfill & GP amplitude (relative flux)\dotfill & $\mathcal{J}(10^{-6},1)$\dotfill &  $%.6f^{+%.6f}_{-%.6f}$ \\'
           %(post['GP_sigma_TESS27'][1],post['GP_sigma_TESS27'][2],post['GP_sigma_TESS27'][3])+"\n")
fout.write(r'~~~~$\rho_{\text{GP,TESS27}}$\dotfill & GP time-scale (days)\dotfill & $\mathcal{J}(10^{-6},10^{3})$\dotfill &  $%.2f^{+%.2f}_{-%.2f}$ \\'
           %(post['GP_rho_TESS27'][1],post['GP_rho_TESS27'][2],post['GP_rho_TESS27'][3])+"\n")


fout.write(r'~~~~$m_{\text{flux,SAAO}}$\dotfill & Offset (relative flux)\dotfill & $\mathcal{N}(0,0.01)$\dotfill &  $%.5f^{+%.5f}_{-%.5f}$ \\'
           %(post['mflux_SAAO'][1],post['mflux_SAAO'][2],post['mflux_SAAO'][3])+"\n")
fout.write(r'~~~~$\sigma_{\text{SAAO}}$\dotfill  & Jitter (ppm)\dotfill & $\mathcal{J}(0.1,10000)$\dotfill & $%.1f^{+%.1f}_{-%.1f}$ \\'
           %(post['sigma_w_SAAO'][1],post['sigma_w_SAAO'][2],post['sigma_w_SAAO'][3])+"\n")

fout.write(r'~~~~$m_{\text{flux,ECAM}}$\dotfill & Offset (relative flux)\dotfill & $\mathcal{N}(0,0.01)$\dotfill &  $%.5f^{+%.5f}_{-%.5f}$ \\'
           %(post['mflux_ECAM'][1],post['mflux_ECAM'][2],post['mflux_ECAM'][3])+"\n")
fout.write(r'~~~~$\sigma_{\text{ECAM}}$\dotfill  & Jitter (ppm)\dotfill & $\mathcal{J}(0.1,10000)$\dotfill & $%.1f^{+%.1f}_{-%.1f}$ \\'
           %(post['sigma_w_ECAM'][1],post['sigma_w_ECAM'][2],post['sigma_w_ECAM'][3])+"\n")

#fout.write(r'~~~~$m_{\text{flux,NGTS}}$\dotfill & Offset (relative flux)\dotfill & $\mathcal{N}(0,0.01)$\dotfill &  $%.5f^{+%.5f}_{-%.5f}$ \\'
#           %(post['mflux_NGTS'][1],post['mflux_NGTS'][2],post['mflux_NGTS'][3])+"\n")
#fout.write(r'~~~~$\sigma_{\text{NGTS}}$\dotfill & Jitter (ppm)\dotfill & $\mathcal{J}(0.1,10000)$\dotfill & $%.2f^{+%.2f}_{-%.2f}$ \\'
#           %(post['sigma_w_NGTS'][1],post['sigma_w_NGTS'][2],post['sigma_w_NGTS'][3])+"\n")

#fout.write(r'~~~~$\sigma_{\text{GP,NGTS}}$\dotfill & GP amplitude (relative flux)\dotfill & $\mathcal{J}(10^{-6},1)$\dotfill &  $%.5f^{+%.5f}_{-%.5f}$ \\'
#           %(post['GP_sigma_NGTS'][1],post['GP_sigma_NGTS'][2],post['GP_sigma_NGTS'][3])+"\n")
#fout.write(r'~~~~$\rho_{\text{GP,NGTS}}$\dotfill & GP time-scale (days)\dotfill & $\mathcal{J}(10^{-6},10^{3})$\dotfill &  $%.3f^{+%.3f}_{-%.3f}$ \\'
#           %(post['GP_rho_NGTS'][1],post['GP_rho_NGTS'][2],post['GP_rho_NGTS'][3])+"\n")


fout.write(r'~~~~$\mu _{\text{CORALIE}}$\dotfill & Systemic RV offset (\kms)\dotfill & $\mathcal{U}(-100,100)$\dotfill & $%.4f^{+%.4f}_{-%.4f}$ \\'
           %(post['mu_CORALIE'][1],post['mu_CORALIE'][2],post['mu_CORALIE'][3])+"\n")
fout.write(r'~~~~$\sigma_{\text{CORALIE}}$\dotfill  & Jitter (\ms)\dotfill &  $\mathcal{J}(0.01,200)$\dotfill & $%.2f^{+%.2f}_{-%.2f}$ \\'
           %(post['sigma_w_CORALIE'][1]*1000.0,post['sigma_w_CORALIE'][2]*1000.0,post['sigma_w_CORALIE'][3]*1000.0)+"\n")

fout.write(r'~~~~$\mu_{\text{HARPS}}$\dotfill & Systemic RV offset (\kms)\dotfill &  $\mathcal{U}(-100,100)$\dotfill & $%.4f^{+%.4f}_{-%.4f}$ \\'
           %(post['mu_HARPS'][1],post['mu_HARPS'][2],post['mu_HARPS'][3])+"\n")
fout.write(r'~~~~$\sigma_{\text{HARPS}}$\dotfill  & Jitter (\ms)\dotfill &  $\mathcal{J}(0.01,200)$\dotfill & $%.2f^{+%.2f}_{-%.2f}$ \\'
           %(post['sigma_w_HARPS'][1]*1000.0,post['sigma_w_HARPS'][2]*1000.0,post['sigma_w_HARPS'][3]*1000.0)+"\n")

fout.write(r''+"\n")
fout.write(r'\smallskip\\\multicolumn{2}{l}{\underline{Modeled Physical Parameters:}}&\smallskip\\'+"\n")
fout.write(r''+"\n")

fout.write(r'~~~~$P$ \dotfill & Period (days) \dotfill & $\mathcal{U}(75.12$\pm$0.2)$\dotfill & $%.5f^{+%.5f}_{-%.5f}$ \\'
           %(post['P_p1'][1],post['P_p1'][2],post['P_p1'][3])+"\n")
fout.write(r'~~~~$T_0$ \dotfill & Time of transit center (BJD$_{\text{TDB}}$) \dotfill & $\mathcal{U}(2458679.3$\pm$0.2)$\dotfill & $%.4f^{+%.4f}_{-%.4f}$ \\'
           %(post['t0_p1'][1]+2400000,post['t0_p1'][2],post['t0_p1'][3])+"\n")
fout.write(r'~~~~$K$\dotfill & Radial velocity semi-amplitude (\ms) \dotfill & $\mathcal{U}(0,1000)$\dotfill & $%.2f^{+%.2f}_{-%.2f}$ \\'
           %(post['K_p1'][1]*1000.0,post['K_p1'][2]*1000.0,post['K_p1'][3]*1000.0)+"\n")
fout.write(r'~~~~$e$\dotfill & Eccentricity of the orbit \dotfill & $\mathcal{B}(0.867,3.03)$\dotfill & $%.3f^{+%.3f}_{-%.3f}$ \\'
           %(post['ecc_p1'][1],post['ecc_p1'][2],post['ecc_p1'][3])+"\n")
fout.write(r'~~~~$\omega$\dotfill & Argument of periastron (deg) \dotfill & $\mathcal{U}(0,360)$\dotfill & $%.1f^{+%.1f}_{-%.1f}$ \\'
           %(post['omega_p1'][1],post['omega_p1'][2],post['omega_p1'][3])+"\n")
fout.write(r'~~~~$r_1$\dotfill & Parametrization for p and b \dotfill & $\mathcal{U}(0,1)$\dotfill & $%.4f^{+%.4f}_{-%.4f}$ \\'
           %(post['r1_p1'][1],post['r1_p1'][2],post['r1_p1'][3])+"\n")
fout.write(r'~~~~$r_2$\dotfill & Parametrization for p and b \dotfill & $\mathcal{U}(0,1)$\dotfill & $%.4f^{+%.4f}_{-%.4f}$ \\'
           %(post['r2_p1'][1],post['r2_p1'][2],post['r2_p1'][3])+"\n")
fout.write(r'~~~~$\rho _{*}$\dotfill & Stellar density (g\,cm$^{-3}$)  \dotfill & $\mathcal{N}(1.07,0.13)$\dotfill & $%.3f^{+%.3f}_{-%.3f}$ \\'
           %(post['rho'][1]/1000.0,post['rho'][2]/1000.0,post['rho'][3]/1000.0)+"\n")

fout.write(r''+"\n")
fout.write(r'\smallskip\\\multicolumn{2}{l}{\underline{Derived Planet Parameters$^{**}$:}}&\smallskip\\'+"\n")
fout.write(r''+"\n")

fout.write(r'~~~~$i$\dotfill & Inclination (deg) \dotfill & \dotfill & $%.3f^{+%.3f}_{-%.3f}$ \\'
           %(post['inc_p1'][1],post['inc_p1'][2],post['inc_p1'][3])+"\n")
fout.write(r'~~~~$p=R_p/R_{\star}$\dotfill & Planet-to-star radius ratio \dotfill & \dotfill & $%.4f^{+%.4f}_{-%.4f}$ \\'
           %(post['p_p1'][1],post['p_p1'][2],post['p_p1'][3])+"\n")
fout.write(r'~~~~$b$\dotfill &Impact parameter of the orbit \dotfill & \dotfill & $%.3f^{+%.3f}_{-%.3f}$ \\'
           %(post['b_p1'][1],post['b_p1'][2],post['b_p1'][3])+"\n")
fout.write(r'~~~~$a$\dotfill & Semi-major axis (AU)  \dotfill & \dotfill & $%.3f^{+%.3f}_{-%.3f}$ \\'
           #%(post['a_p1'][1],post['a_p1'][2],post['a_p1'][3])+"\n")
           %(a_AU_value,a_AU_value_up,a_AU_value_low)+"\n")
fout.write(r'~~~~$M_p$\dotfill & Planetary mass (\mjup)  \dotfill & \dotfill & $%.2f^{+%.2f}_{-%.2f}$ \\'
           %(Mp_value,Mp_value_up,Mp_value_low)+"\n")
fout.write(r'~~~~$R_p$\dotfill & Planetary radius (\rjup)  \dotfill & \dotfill & $%.3f^{+%.3f}_{-%.3f}$ \\'
           %(Rp_value,Rp_value_up,Rp_value_low)+"\n")
fout.write(r'~~~~$\rho _p$\dotfill & Planetary density (g\,cm$^{-3}$)  \dotfill & \dotfill & $%.2f^{+%.2f}_{-%.2f}$ \\'
           %(planet_density_value,planet_density_value_up,planet_density_value_low)+"\n")
fout.write(r'~~~~$S$\dotfill & Insolation ($S_{\oplus}$)  \dotfill & \dotfill & $%.1f^{+%.1f}_{-%.1f}$ \\'
           %(S_value,S_value_up,S_value_low)+"\n")
fout.write(r'~~~~$T_{eq}$\dotfill & Equilibrium Temperature ($K$) \dotfill & \dotfill & $%i^{+%i}_{-%i}$ \\'
           %(Teq_value,Teq_value_up,Teq_value_low)+"\n")

fout.write(r''+"\n")
fout.write(r'\end{table*}'+ "\n")
fout.close()