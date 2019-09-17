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


#options = ({'DoLateRadTruncation':False})


#GET RESULTS FOR SEVERAL CASES
pars = camb.CAMBparams()
#pars.set_accuracy(lAccuracyBoost = 10)
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency

pars.set_cosmology(H0=67.0, ombh2=0.02237, omch2= 0.12011, mnu=0.06, omk=0, tau=0.0543,
num_bins=1, smooth_factor=10, void_model=2, zbins0=3000.0, qbins0 = 0.0)

pars.InitPower.set_params(ns=0.965, r=0, As=2e-9)
pars.set_for_lmax(2500, lens_potential_accuracy=0);

results_dhost = camb.get_results(pars)
pars.NonLinear = model.NonLinear_both

results_dhost.calc_power_spectra(pars)

powers =results_dhost.get_cmb_power_spectra(pars, CMB_unit='muK')

z  = np.linspace(0,3000,10000)

O_m = [results_dhost.rhoc_of_z(zi)/results_dhost.h_of_z(zi)**2./3. for zi in z]
O_v = [results_dhost.rhov_of_z(zi)/results_dhost.h_of_z(zi)**2./3. for zi in z]

np.savetxt('initsolvertest_om_lcdm.dat', O_m)
np.savetxt('initsolvertest_ov_lcdm.dat', O_v)

