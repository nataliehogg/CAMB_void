##matplotlib inline
import sys, platform, os
from matplotlib import pyplot as plt
import numpy as np
print('Using CAMB installed at %s'%(os.path.realpath(os.path.join(os.getcwd(),'.'))))
#uncomment this if you are running remotely and want to keep in synch with repo changes
#if platform.system()!='Windows':
#    !cd $HOME/git/camb; git pull github master; git log -1
sys.path.insert(0,os.path.realpath(os.path.join(os.getcwd(),'..')))
import camb
from camb import model, initialpower
print('CAMB version: %s '%camb.__version__)


#GET RESULTS FOR SEVERAL CASES

pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=70.0, ombh2=0.0226, omch2=0.112, mnu=0.06, omk=0, tau=0.06, void_interaction=7, rhov_t=1, num_bins=1, smooth_factor=10,
void_model =2, zbins0=1.0, qbins0=0.0)
pars.InitPower.set_params(ns=0.965, r=0, As=2e-9)
pars.set_for_lmax(2500, lens_potential_accuracy=0);
results_dhost = camb.get_results(pars)

pars.NonLinear = model.NonLinear_both
results_dhost.calc_power_spectra(pars)
kh_nonlin_LCDM, z_nonlin_LCDM, pk_nonlin_LCDM = results_dhost.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)

z  = np.linspace(0,4,100)

O_m_LCDM = [results_dhost.rhoc_of_z(zi)/results_dhost.h_of_z(zi)**2./3. for zi in z]
O_v_LCDM = [results_dhost.rhov_of_z(zi)/results_dhost.h_of_z(zi)**2./3. for zi in z]
H_LCDM   = results_dhost.hubble_parameter(z)


#STARTING POINT
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=70.0, ombh2=0.0226, omch2=0.112, mnu=0.06, omk=0, tau=0.06, void_interaction=7, rhov_t=1, num_bins=1, smooth_factor=10,
void_model =2, zbins0=1.0, qbins0=0.1)
pars.InitPower.set_params(ns=0.965, r=0, As=2e-9)
pars.set_for_lmax(2500, lens_potential_accuracy=0);
results_dhost = camb.get_results(pars)

pars.NonLinear = model.NonLinear_both
results_dhost.calc_power_spectra(pars)
kh_nonlin, z_nonlin, pk_nonlin = results_dhost.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)

z  = np.linspace(0,4,100)

O_m = [results_dhost.rhoc_of_z(zi)/results_dhost.h_of_z(zi)**2./3. for zi in z]
O_v = [results_dhost.rhov_of_z(zi)/results_dhost.h_of_z(zi)**2./3. for zi in z]
H   = results_dhost.hubble_parameter(z)

pars.set_cosmology(H0=70.0, ombh2=0.0226, omch2=0.112, mnu=0.06, omk=0, tau=0.06, void_interaction=7, rhov_t=1, num_bins=1, smooth_factor=10,
void_model =2, zbins0=1.0, qbins0=-0.1)
pars.InitPower.set_params(ns=0.965, r=0, As=2e-9)
pars.set_for_lmax(2500, lens_potential_accuracy=0);
results_dhost = camb.get_results(pars)

pars.NonLinear = model.NonLinear_both
results_dhost.calc_power_spectra(pars)
kh_nonlin_2, z_nonlin_2, pk_nonlin_2 = results_dhost.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)

z  = np.linspace(0,4,100)

O_m_2 = [results_dhost.rhoc_of_z(zi)/results_dhost.h_of_z(zi)**2./3. for zi in z]
O_v_2 = [results_dhost.rhov_of_z(zi)/results_dhost.h_of_z(zi)**2./3. for zi in z]
H_2   = results_dhost.hubble_parameter(z)

#print(rho_m, rho_v)

#Plots and stuff
plt.figure()
plt.subplot(111)
plt.plot(z, O_v_LCDM, color='black', label=r'$\Lambda$CDM')
plt.plot(z, O_v, color='#FFB300', label=r'$q_V =0.1$')
plt.plot(z, O_v_2, color='#8E001C', label=r'$q_V =-0.1$')
plt.plot(z, O_m_LCDM, color='black',ls='--')
plt.plot(z, O_m, color='#FFB300',ls='--')
plt.plot(z, O_m_2, color='#8E001C',ls='--')
plt.ylabel(r'$\Omega(z)$', fontsize=12)
plt.xlabel(r'$z$', fontsize=12)
plt.xlim(0,3)
plt.ylim(0,1)
plt.legend(loc='center right')#, fontsize='small');
plt.savefig('omegas.pdf')


#Plots and stuff
plt.figure()
plt.subplot(111)
plt.loglog(kh_nonlin_LCDM, pk_nonlin_LCDM[0,:], color='black', label=r'$\Lambda$CDM')
plt.loglog(kh_nonlin, pk_nonlin[0,:], color='#FFB300', label=r'$q_V =0.1$')
plt.loglog(kh_nonlin_2, pk_nonlin_2[0,:], color='#8E001C', label=r'$q_V =-0.1$')
plt.ylabel(r'$P(k,z=0)$', fontsize=12)
plt.xlabel(r'$k/h$', fontsize=12)
plt.legend(loc='upper right')#, fontsize='small');
plt.savefig('matterpower.pdf')

