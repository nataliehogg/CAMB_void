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


#STARTING POINT
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=70.0, ombh2=0.0226, omch2=0.112, mnu=0.06, omk=0, tau=0.06, void_interaction=7, rhov_t=1, num_bins=1, smooth_factor=10,
void_model =2, zbins1=1.0, qbins1=1)
pars.InitPower.set_params(ns=0.965, r=0, As=2e-9)
pars.set_for_lmax(2500, lens_potential_accuracy=0);
#pars.minimizeme = False
#calculate results for these parameters
results_dhost = camb.get_results(pars)

z  = np.linspace(0,4,100)

O_m = [results_dhost.rhoc_of_z(zi)/results_dhost.h_of_z(zi)**2./3. for zi in z]
O_v = [results_dhost.rhov_of_z(zi)/results_dhost.h_of_z(zi)**2./3. for zi in z]

#print(rho_m, rho_v)

#Plots and stuff
plt.plot(z, O_m, color='black', label=r'$\Omega_m$')
plt.plot(z, O_v, color='#8E001C', label=r'$\Omega_v$')
plt.ylabel(r'$\Omega$')
plt.xlim(0,4)
plt.ylim(0,1)
# plt.yscale('log')
plt.legend(loc='best', fontsize='small');
plt.savefig('newinteractiondensities.pdf')
#plt.show()
