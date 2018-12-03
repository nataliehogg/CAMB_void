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
pars.set_cosmology(H0=67.32, ombh2=0.02237, omch2= 0.12011, mnu=0.06, omk=0, tau=0.0543, void_interaction=7, rhov_t=1, num_bins=1, smooth_factor=10,
void_model =2, zbins0=1.0, qbins0=0.0)
pars.InitPower.set_params(ns=0.96605, r=0, As=2.0989e-9)
pars.set_for_lmax(2500, lens_potential_accuracy=0);
results_dhost = camb.get_results(pars)
powers =results_dhost.get_cmb_power_spectra(pars, CMB_unit='muK')
pars.NonLinear = model.NonLinear_both
results_dhost.calc_power_spectra(pars)
kh_nonlin_LCDM, z_nonlin_LCDM, pk_nonlin_LCDM = results_dhost.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)

z  = np.linspace(0,4,100)

O_m_LCDM = [results_dhost.rhoc_of_z(zi)/results_dhost.h_of_z(zi)**2./3. for zi in z]
O_v_LCDM = [results_dhost.rhov_of_z(zi)/results_dhost.h_of_z(zi)**2./3. for zi in z]
H_LCDM   = results_dhost.hubble_parameter(z)
DL_LCDM = results_dhost.luminosity_distance(z)
totCL_LCDM=powers['total']
ls_LCDM = np.arange(totCL_LCDM.shape[0])
print('CosmoMC theta_MC parameter: %s'%results_dhost.cosmomc_theta())

#STARTING POINT
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=55, ombh2=0.02237, omch2= 0.20386, mnu=0.06, omk=0, tau=0.0543, void_interaction=7, rhov_t=1, num_bins=1, smooth_factor=10,
void_model =2, zbins0=3000.0, qbins0=1.6)
pars.InitPower.set_params(ns=0.96605, r=0, As=2.0989e-9)
pars.set_for_lmax(2500, lens_potential_accuracy=0);
results_dhost = camb.get_results(pars)

pars.NonLinear = model.NonLinear_both
results_dhost.calc_power_spectra(pars)
powers =results_dhost.get_cmb_power_spectra(pars, CMB_unit='muK')
kh_nonlin, z_nonlin, pk_nonlin = results_dhost.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
z  = np.linspace(0,4,100)

O_m = [results_dhost.rhoc_of_z(zi)/results_dhost.h_of_z(zi)**2./3. for zi in z]
O_v = [results_dhost.rhov_of_z(zi)/results_dhost.h_of_z(zi)**2./3. for zi in z]
H   = results_dhost.hubble_parameter(z)
DL = results_dhost.luminosity_distance(z)
totCL=powers['total']
ls = np.arange(totCL.shape[0])
print('CosmoMC theta_MC parameter: %s'%results_dhost.cosmomc_theta())

pars.set_cosmology(H0=68, ombh2=0.02237, omch2= 0.101833, mnu=0.06, omk=0, tau=0.0543, void_interaction=7, rhov_t=1, num_bins=1, smooth_factor=10,
void_model =2, zbins0=3000.0, qbins0=-0.1)
pars.InitPower.set_params(ns=0.96605, r=0, As=2.0989e-9)
pars.set_for_lmax(2500, lens_potential_accuracy=0);
results_dhost = camb.get_results(pars)
powers =results_dhost.get_cmb_power_spectra(pars, CMB_unit='muK')
pars.NonLinear = model.NonLinear_both
results_dhost.calc_power_spectra(pars)
kh_nonlin_2, z_nonlin_2, pk_nonlin_2 = results_dhost.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)

z  = np.linspace(0,4,100)

O_m_2 = [results_dhost.rhoc_of_z(zi)/results_dhost.h_of_z(zi)**2./3. for zi in z]
O_v_2 = [results_dhost.rhov_of_z(zi)/results_dhost.h_of_z(zi)**2./3. for zi in z]
H_2   = results_dhost.hubble_parameter(z)
DL_2 = results_dhost.luminosity_distance(z)
totCL_2=powers['total']
ls_2= np.arange(totCL_2.shape[0])

#print(rho_m, rho_v)
print('CosmoMC theta_MC parameter: %s'%results_dhost.cosmomc_theta())
print(ls_2)

plt.figure()
plt.subplot(111)
plt.plot(z, H_LCDM, color='black', label=r'$\Lambda$CDM, $H_0=67.36$')
plt.plot(z, H, color='#FFB300', label=r'$q_V =0.1,\ H_0=71$')
plt.plot(z, H_2, color='#8E001C', label=r'$q_V =-0.1,\ H_0=65$')
plt.ylabel(r'$H(z)$', fontsize=12)
plt.xlabel(r'$z$', fontsize=12)
#plt.xlim(0,3)
#plt.ylim(0,1)
plt.legend(loc='lower right')#, fontsize='small');
#plt.savefig('hubble.pdf')
#plt.show()

data = np.loadtxt('data.dat')

plt.figure()
plt.subplot(111)
plt.plot(ls[2:], totCL_LCDM[2:,0], color='black', label=r'$\Lambda$CDM, $H_0=67.36$')
plt.plot(ls[2:], totCL[2:,0], color='#FFB300', label=r'Planck')
plt.plot(ls[2:], totCL_2[2:,0], color='#8E001C', label=r'Planck+LSS')
plt.errorbar(data[:,0], data[:,2], yerr=[data[:,3],data[:,4]], marker='*', markersize=2, ecolor='black',ls='')
plt.ylabel(r'$\ell(\ell+1)/2\pi C_\ell^{TT}$', fontsize=12)
plt.xlabel(r'$\ell$', fontsize=12)
plt.xscale('log')
#plt.xlim(0,3)
#plt.ylim(0,1)
plt.legend(loc='upper left')#, fontsize='small');
plt.savefig('TTspectra_bestfits.pdf')
#plt.show()



plt.figure()
plt.subplot(111)
plt.loglog(kh_nonlin_LCDM, pk_nonlin_LCDM[0,:], color='black', label=r'$\Lambda$CDM')
plt.loglog(kh_nonlin, pk_nonlin[0,:], color='#FFB300', label=r'Planck')
plt.loglog(kh_nonlin_2, pk_nonlin_2[0,:], color='#8E001C', label=r'Planck+LSS')
plt.ylabel(r'$P(k,z=0)$', fontsize=12)
plt.xlabel(r'$k/h$', fontsize=12)
plt.legend(loc='upper right')#, fontsize='small');
plt.savefig('matterpower_bestfits.pdf')
#plt.show()
