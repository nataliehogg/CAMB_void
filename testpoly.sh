#!/bin/env bash
##comment out lines by adding at least two `#' at the beginning
#SBATCH --job-name=testpoly
#SBATCH --account=silvestri
#SBATCH --partition=computation
#SBATCH --output=/marisdata/martinelli/MY_PROJECTS/ONGOING/void/polyversion/CAMB_void/logfiles/%x.out
#SBATCH --error=/marisdata/martinelli/MY_PROJECTS/ONGOING/void/polyversion/CAMB_void/logfiles/%x.err
#SBATCH --time=1-00:00:00
#SBATCH --mem=48442
#SBATCH --ntasks=4
# run hostname for instance

ulimit -s unlimited

mpirun -n 4 ./cosmomc parameter_files/Planck_16binsCP.ini

