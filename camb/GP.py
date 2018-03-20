from sklearn.gaussian_process import GaussianProcessRegressor
import sys
from sklearn.gaussian_process.kernels import RBF
import numpy as np
import matplotlib.pyplot as plt
import argparse



parser = argparse.ArgumentParser(description='Input parameters from CAMB')
parser.add_argument('--inired',metavar='zini',  type=float, nargs='+',
                   help='initial redshift')
parser.add_argument('--endred',metavar='zend',  type=float, nargs='+',
                   help='end redshift')
parser.add_argument('--ODEsteps',metavar='ODEsteps',  type=int, nargs='+',
                   help='number of steps for the ODE solver')
parser.add_argument('--redshifts',metavar='z',  type=float, nargs='*',default=[],
                   help='values of redshifts')
parser.add_argument('--couplings',metavar='q',  type=float, nargs='*',default=[],
                   help='values of couplings')
parser.add_argument('--l',metavar='l',  type=float, nargs='+',
                   help='correlation length')
parser.add_argument('--outfile', nargs='+', type=str ,default=sys.stdout)
args = parser.parse_args()


#Training points
inired = args.inired[0]
endred = args.endred[0]
ODEsteps = args.ODEsteps[0]
q = np.array(args.couplings)
l = args.l[0]
filename = args.outfile[0]

z = np.array(args.redshifts)#z_edge[:-1] + np.diff(z_edge)/2 #NH z is now the redshift in the middle of each bin

# Generation of the Gaussian Process
gp = GaussianProcessRegressor(kernel=RBF(l, (l, l)))

#Fit --> Training
g = gp.fit(z[:, np.newaxis], q)

#Plotting points (if log use np.logspace)
z_sampling = np.linspace(inired, endred, ODEsteps)

#Predict points
q_pred, sigma = gp.predict(z_sampling[:, np.newaxis], return_std=True)

np.savetxt(filename, np.array([z_sampling, q_pred]).T, fmt="%15.8e")


exit()
