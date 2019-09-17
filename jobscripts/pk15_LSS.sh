#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=20:00
#SBATCH --job-name=plik18_fixfid
#SBATCH --output=jobscripts/plik18_fixfid.out
#SBATCH --error=jobscripts/plik18_fixfid.err
#SBATCH -p sciama4.q
#SBATCH --mail-type=ALL
#SBATCH --mail-user=natalie.hogg@port.ac.uk
#SBATCH -D /users/hoggn/building_17bin_plik18

# Configure modules and load modules
module load /opt/apps/etc/modules/core/services/s3cmd
module load /opt/apps/etc/modules/core/system/intel64
module load /opt/apps/etc/modules/compilers/gnu_comp/5.4.0
module load /opt/apps/etc/modules/compilers/intel_comp/2019.2
module load /opt/apps/etc/modules/libraries/cfitsio/3.47
module load /opt/apps/etc/modules/applications/anaconda3/2019.03
module load /opt/apps/etc/modules/mpi/openmpi/4.0.1

mpirun -np 8 ./cosmomc inifiles/plik18_fixfid.ini > chains/plik18_fixfid.log


