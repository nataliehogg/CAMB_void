#!/bin/bash
#PBS -l nodes=8:ppn=2
#PBS -l walltime=48:00:00
#PBS -N 4bin_Planck_model1
#PBS -o 4bin_Planck_model1
#PBS -e 4bin_Planck_model1
#PBS -m abe
#PBS -M natalie.hogg@port.ac.uk
#PBS -V
#PBS -q cluster.q
#PBS -j oe
#PBS -d /users/hoggn/CosmoMC_void

# Configure modules and load modules
module load compilers/intel/intel-64.v15.2.164
module load mpi/openmpi/1.8.4/intel-64.v15.2.164
module load libs/mkl/2015.2.164/intel-64.v15.2.164
module load apps/python/2.7.9/intel-64.v15.2.164

mpirun -np 4 ./cosmomc parameter_files/4bin_Planck_model1.ini
