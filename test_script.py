##matplotlib inline
import sys, platform, os
from matplotlib import pyplot as plt
import numpy as np
print('Using CAMB installed at %s'%(os.path.realpath(os.path.join(os.getcwd(),'./'))))
#uncomment this if you are running remotely and want to keep in synch with repo changes
#if platform.system()!='Windows':
#    !cd $HOME/git/camb; git pull github master; git log -1
sys.path.insert(0,os.path.realpath(os.path.join(os.getcwd(),'./')))
import camb
from camb import model, initialpower
print('CAMB version: %s '%camb.__version__)


#GET RESULTS FOR SEVERAL CASES


#MINIMIZED DHOST
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=70.0, ombh2=0.0226, omch2=0.112, mnu=0.06, omk=0, tau=0.06, num_bins = 3, ending_z = 2.0, zbins = [1.0,1.5,2.0], qbins = [0.,0.,0.])
pars.InitPower.set_params(ns=0.965, r=0, As=2e-9)
pars.set_for_lmax(2500, lens_potential_accuracy=0);
#calculate results for these parameters
results_dhost = camb.get_results(pars)

z  = np.linspace(0,4,100)
DA_dhost = results_dhost.angular_diameter_distance(z)
H_dhost  = results_dhost.hubble_parameter(z)
DL_dhost = results_dhost.luminosity_distance(z)


##LCDM LIMIT (setting inired=0)
#Set up a new set of parameters for CAMB
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
pars.set_cosmology(H0=70.0, ombh2=0.0226, omch2=0.112, mnu=0.06, omk=0, tau=0.06, void_model = 0, num_bins = 2, ending_z = 2.0, zbins = [1.0,2.0], qbins = [0.0,0.0])
pars.InitPower.set_params(ns=0.965, r=0, As=2e-9)
pars.set_for_lmax(2500, lens_potential_accuracy=0);
pars.minimizeme = False
#calculate results for these parameters
results_lcdm = camb.get_results(pars)


DA = results_lcdm.angular_diameter_distance(z)
H  = results_lcdm.hubble_parameter(z)
DL = results_lcdm.luminosity_distance(z)
#----------------------------------------------------------



#Plots and stuff
plt.plot(z, DA, color='#8E001C', label=r'$\Lambda$CDM')
plt.plot(z, DA_dhost, color='#FFB300', label=r'test void')
plt.xlabel('$z$')
plt.ylabel(r'$D_A /\rm{Mpc}$')
plt.title('Angular diameter distance')
plt.xlim([0,4])
plt.ylim([0,2500]);
plt.legend(loc='lower right', fontsize='small');
plt.savefig('da_dhost.pdf')
plt.show()


#Plots and stuff
plt.plot(z, H, color='#8E001C', label=r'$\Lambda$CDM')
plt.plot(z, H_dhost, color='#FFB300', label=r'test void')
plt.xlabel('$z$')
plt.ylabel(r'$H(z)\ {\rm km}/{\rm s}/{\rm Mpc}$')
plt.title('Hubble parameter')
#plt.xlim([0,2500]);
plt.legend(loc='lower right', fontsize='small');
plt.savefig('H_dhost.pdf')
plt.show()

#Plots and stuff
plt.plot(z, DL, color='#8E001C', label=r'$\Lambda$CDM')
plt.plot(z, DL_dhost, color='#FFB300', label=r'test void')
plt.xlabel('$z$')
plt.ylabel(r'$D_L /\rm{Mpc}$')
plt.title('Luminosity distance')
plt.xlim([0,4])
#plt.ylim([0,2500]);
plt.legend(loc='lower right', fontsize='small');
plt.savefig('lumdis_dhost.pdf')
plt.show()
