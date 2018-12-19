##matplotlib inline
import sys, platform, os
from matplotlib import pyplot as plt
import numpy as np
print('Using CAMB installed at %s'%(os.path.realpath(os.path.join(os.getcwd(),'./pycamb'))))
#uncomment this if you are running remotely and want to keep in synch with repo changes
#if platform.system()!='Windows':
#    !cd $HOME/git/camb; git pull github master; git log -1
sys.path.insert(0,os.path.realpath(os.path.join(os.getcwd(),'./pycamb')))
import camb
from camb import model, initialpower
print('CAMB version: %s '%camb.__version__)


#GET RESULTS FOR SEVERAL CASES



#THIS IS LCDM WITH PLANCK BEST FIT
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=67.32, ombh2=0.02237, omch2= 0.12011, mnu=0.06, omk=0, tau=0.0543, void_interaction=7, rhov_t=1, num_bins=1, smooth_factor=10,
void_model =2, zbins0=10., qbins0=0.)
pars.InitPower.set_params(ns=0.96605, r=0, As=2.0989e-9)
pars.set_for_lmax(2500, lens_potential_accuracy=0);
results_dhost = camb.get_results(pars)
powers =results_dhost.get_cmb_power_spectra(pars, CMB_unit='muK')
pars.NonLinear = model.NonLinear_both
results_dhost.calc_power_spectra(pars)
kh_nonlin_LCDM, z_nonlin_LCDM, pk_nonlin_LCDM = results_dhost.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)

z  = np.linspace(0,4,100)

rhoc_LCDM = [results_dhost.rhoc_of_z(zi) for zi in z]
rhov_LCDM = [results_dhost.rhov_of_z(zi) for zi in z]
O_m_LCDM = [results_dhost.rhoc_of_z(zi)/results_dhost.h_of_z(zi)**2./3. for zi in z]
O_v_LCDM = [results_dhost.rhov_of_z(zi)/results_dhost.h_of_z(zi)**2./3. for zi in z]
H_LCDM   = results_dhost.hubble_parameter(z)
totCL_LCDM=powers['total']
ls_LCDM = np.arange(totCL_LCDM.shape[0])
print('CosmoMC theta_MC parameter: %s'%results_dhost.cosmomc_theta())


#THIS IS THE VVE BEST FIT
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=6.5250950E+01, ombh2=2.2374700E-02, omch2=2.0481860E-01, mnu=0.06, omk=0, tau=1.8052370E-02, void_interaction=8, rhov_t=1, num_bins=1, smooth_factor=10,
void_model =2, zbins0=1.0450980E-01, qbins0=-2.1193660)
pars.InitPower.set_params(ns=9.6729770E-01, r=0, As=1.e-10*np.exp(2.9739370))#2.0989e-9)
pars.set_for_lmax(2500, lens_potential_accuracy=0);
results_dhost = camb.get_results(pars)
powers =results_dhost.get_cmb_power_spectra(pars, CMB_unit='muK')
pars.NonLinear = model.NonLinear_both
results_dhost.calc_power_spectra(pars)
kh_nonlin, z_nonlin, pk_nonlin = results_dhost.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)

z  = np.linspace(0,4,100)

rhoc = [results_dhost.rhoc_of_z(zi) for zi in z]
rhov = [results_dhost.rhov_of_z(zi) for zi in z]
O_m  = [results_dhost.rhoc_of_z(zi)/results_dhost.h_of_z(zi)**2./3. for zi in z]
O_v  = [results_dhost.rhov_of_z(zi)/results_dhost.h_of_z(zi)**2./3. for zi in z]
H    = results_dhost.hubble_parameter(z)
totCL=powers['total']
ls   = np.arange(totCL.shape[0])
print('CosmoMC theta_MC parameter: %s'%results_dhost.cosmomc_theta())



#THIS IS THE CVAR BEST FIT
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=6.9574850E+01, ombh2=2.2271410E-02, omch2=9.4193890E-02, mnu=0.06, omk=0, tau=7.7410870E-02, void_interaction=7, rhov_t=1, num_bins=1, smooth_factor=10,
void_model =2, zbins0=5.2004670E+00, qbins0=-2.2471370E-01)
pars.InitPower.set_params(ns=9.6667080E-01, r=0, As=1.e-10*np.exp(3.0872150E+00))
pars.set_for_lmax(2500, lens_potential_accuracy=0);
results_dhost = camb.get_results(pars)
powers =results_dhost.get_cmb_power_spectra(pars, CMB_unit='muK')
pars.NonLinear = model.NonLinear_both
results_dhost.calc_power_spectra(pars)
kh_nonlin_2, z_nonlin_2, pk_nonlin_2 = results_dhost.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)

z  = np.linspace(0,4,100)

rhoc_2 = [results_dhost.rhoc_of_z(zi) for zi in z]
rhov_2 = [results_dhost.rhov_of_z(zi) for zi in z]
O_m_2  = [results_dhost.rhoc_of_z(zi)/results_dhost.h_of_z(zi)**2./3. for zi in z]
O_v_2  = [results_dhost.rhov_of_z(zi)/results_dhost.h_of_z(zi)**2./3. for zi in z]
H_2    = results_dhost.hubble_parameter(z)
totCL_2=powers['total']
ls_2   = np.arange(totCL_2.shape[0])
print('CosmoMC theta_MC parameter: %s'%results_dhost.cosmomc_theta())





plt.figure()
plt.subplot(111)
plt.plot(z, rhoc     , color='#FFB300', ls='-' , label=r'VVE best fit')
plt.plot(z, rhov     , color='#FFB300', ls='--')
plt.plot(z, rhoc_2   , color='#8E001C', ls='-' , label=r'Cvar best fit')
plt.plot(z, rhov_2   , color='#8E001C', ls='--')
plt.plot(z, rhoc_LCDM, color='black'  , ls='-' , label=r'$\Lambda$CDM')
plt.plot(z, rhov_LCDM, color='black'  , ls='--')
plt.ylabel(r'$\rho(z)$', fontsize=12)
plt.xlabel(r'$z$', fontsize=12)
plt.xlim(0,3)
#plt.ylim(0,1)
plt.legend(loc='upper left')#, fontsize='small');
#plt.savefig('rhos_postfix.pdf')
plt.show()

plt.figure()
plt.subplot(111)
plt.plot(z, O_m     , color='#FFB300', ls='-' , label=r'VVE best fit')
plt.plot(z, O_v     , color='#FFB300', ls='--')
plt.plot(z, O_m_2   , color='#8E001C', ls='-' , label=r'Cvar best fit')
plt.plot(z, O_v_2   , color='#8E001C', ls='--')
plt.plot(z, O_m_LCDM, color='black'  , ls='-' , label=r'$\Lambda$CDM')
plt.plot(z, O_v_LCDM, color='black'  , ls='--')
plt.ylabel(r'$\rho(z)$', fontsize=12)
plt.xlabel(r'$z$', fontsize=12)
plt.xlim(0,3)
#plt.ylim(0,1)
plt.legend(loc='upper left')#, fontsize='small');
plt.savefig('omega_bestfit.pdf')
plt.show()


plt.figure()
plt.subplot(111)
plt.plot(z, H_LCDM, color='black', label=r'$\Lambda$CDM')
plt.plot(z, H, color='#FFB300', label=r'VVE best fit')
plt.plot(z, H_2, color='#8E001C', label=r'Cvar best fit')
plt.ylabel(r'$H(z)$', fontsize=12)
plt.xlabel(r'$z$', fontsize=12)
#plt.xlim(0,3)
#plt.ylim(0,1)
plt.legend(loc='lower right')#, fontsize='small');
plt.savefig('hubble_bestfit.pdf')
plt.show()

data = np.loadtxt('data.dat')

plt.figure()
plt.subplot(111)
plt.plot(ls[2:], totCL_LCDM[2:,0], color='black', label=r'$\Lambda$CDM')
plt.plot(ls[2:], totCL[2:,0], color='#FFB300', label=r'VVE best fit')
plt.plot(ls_2[2:], totCL_2[2:,0], color='#8E001C', label=r'Cvar best fit')
plt.errorbar(data[:,0], data[:,2], yerr=[data[:,3],data[:,4]], marker='*', markersize=2, ecolor='black',ls='')
plt.ylabel(r'$\ell(\ell+1)/2\pi C_\ell^{TT}$', fontsize=12)
plt.xlabel(r'$\ell$', fontsize=12)
plt.xscale('log')
#plt.xlim(0,3)
#plt.ylim(0,1)
plt.legend(loc='upper left')#, fontsize='small');
plt.savefig('TT_bestfit.pdf')
plt.show()
