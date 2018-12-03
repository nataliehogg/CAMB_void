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
pars.set_cosmology(H0=67.32, ombh2=0.02237, omch2= 0.12011, mnu=0.06, omk=0, tau=0.0543, void_interaction=8, rhov_t=1, num_bins=1, smooth_factor=10,
void_model =2, zbins0=2.0, qbins0=-0.5)
pars.InitPower.set_params(ns=0.96605, r=0, As=2.0989e-9)
pars.set_for_lmax(2500, lens_potential_accuracy=0);
results_dhost = camb.get_results(pars)
powers =results_dhost.get_cmb_power_spectra(pars, CMB_unit='muK')
pars.NonLinear = model.NonLinear_both
results_dhost.calc_power_spectra(pars)
kh_nonlin_LCDM, z_nonlin_LCDM, pk_nonlin_LCDM = results_dhost.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)

z  = np.linspace(0,4,100)

rhoc = [results_dhost.rhoc_of_z(zi) for zi in z]
rhov = [results_dhost.rhov_of_z(zi) for zi in z]
print('CosmoMC theta_MC parameter: %s'%results_dhost.cosmomc_theta())


plt.figure()
plt.subplot(111)
plt.plot(z, rhoc, color='#FFB300', label=r'$\rho_c(z)$')
plt.plot(z, rhov, color='#8E001C', label=r'$V(z)$')
plt.ylabel(r'$\rho(z)$', fontsize=12)
plt.xlabel(r'$z$', fontsize=12)
plt.xlim(0,3)
#plt.ylim(0,1)
plt.legend(loc='upper left')#, fontsize='small');
plt.savefig('rhos_postfix.pdf')
plt.show()

