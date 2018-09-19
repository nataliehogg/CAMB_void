##matplotlib inline
import sys, platform, os
from matplotlib import pyplot as plt
import numpy as np
print('Using CAMB installed at %s'%(os.path.realpath(os.path.join(os.getcwd(),'../../CAMB_original'))))
#uncomment this if you are running remotely and want to keep in synch with repo changes
#if platform.system()!='Windows':
#    !cd $HOME/git/camb; git pull github master; git log -1
sys.path.insert(0,os.path.realpath(os.path.join(os.getcwd(),'..')))
import camb
from camb import model, initialpower
print('CAMB version: %s '%camb.__version__)


#MINIMIZED VOID
pars = camb.CAMBparams()
#This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
#pars.set_cosmology(H0=70.0, ombh2=0.0226, omch2=0.112, mnu=0.06, omk=0, tau=0.06)
pars.set_cosmology(H0=70.0, ombh2=0.0226, omch2=0.112, mnu=0.06, omk=0, tau=0.06, void_model=2, num_bins = 4,
        smooth_factor = 10, bin_redshift_1 = 0.3,bin_q_1 = -0.1,bin_redshift_2 = 0.9,bin_q_2 = 0.3,bin_redshift_3 = 2.5,bin_q_3 = -1.4,bin_redshift_4 = 10,bin_q_4 = 0.0, correlation_length = 0.5,
        ending_z   = 10, ODEsteps   = 10000)
pars.InitPower.set_params(ns=0.965, r=0, As=2e-9)
pars.set_for_lmax(2500, lens_potential_accuracy=0);
#pars.minimizeme = False
#calculate results for these parameters
results_VOID = camb.get_results(pars)

z  = np.linspace(0,4,100)
DA_VOID = results_VOID.angular_diameter_distance(z)
H_VOID  = results_VOID.hubble_parameter(z)
DL_VOID = results_VOID.luminosity_distance(z)

#----------------------------------------------------------



#Plots and stuff
plt.plot(z, DA_VOID, color='#FFB300', label=r'VOID ($z_{\rm ini}=10$)')
plt.xlabel('$z$')
plt.ylabel(r'$D_A /\rm{Mpc}$')
plt.title('Angular diameter distance')
plt.xlim([0,4])
plt.ylim([0,2500]);
plt.legend(loc='lower right', fontsize='small');
plt.savefig('da_VOID.pdf')
plt.show()


#Plots and stuff
plt.plot(z, H_VOID, color='#FFB300', label=r'VOID ($z_{\rm ini}=10$)')
plt.xlabel('$z$')
plt.ylabel(r'$H(z)\ {\rm km}/{\rm s}/{\rm Mpc}$')
plt.title('Hubble parameter')
#plt.xlim([0,2500]);
plt.legend(loc='lower right', fontsize='small');
plt.savefig('H_VOID.pdf')
plt.show()

#Plots and stuff
plt.plot(z, DL_VOID, color='#FFB300', label=r'VOID ($z_{\rm ini}=10$)')
plt.xlabel('$z$')
plt.ylabel(r'$D_L /\rm{Mpc}$')
plt.title('Luminosity distance')
plt.xlim([0,4])
#plt.ylim([0,2500]);
plt.legend(loc='lower right', fontsize='small');
plt.savefig('lumdis_VOID.pdf')
plt.show()
