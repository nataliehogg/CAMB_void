#!/bin/env bash
##comment out lines by adding at least two `#' at the beginning
#SBATCH --job-name=getdist
#SBATCH --account=silvestri
#SBATCH --partition=computation
#SBATCH --output=/marisdata/martinelli/MY_PROJECTS/ONGOING/void/newcosmo/cosmomc/%x.out
#SBATCH --error=/marisdata/martinelli/MY_PROJECTS/ONGOING/void/newcosmo/cosmomc/%x.err
#SBATCH --time=1-00:00:00
#SBATCH --mem=16681
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
export OMP_NUM_THREADS=5
export MKL_NUM_THREADS=5

#wc chains/*.txt

python plotters/plotter_triangle_Cfix.py
python plotters/plotter_triangle_Cvar.py  
python plotters/plotter_triangle_VVE.py
python plotters/plotter_triangle_VVE_standard.py


python python/GetDist.py chains/pk15_Cfix
python python/GetDist.py chains/pk15_LSS_Cfix
python python/GetDist.py chains/pk15_Cvar
python python/GetDist.py chains/pk15_LSS_Cvar
python python/GetDist.py chains/pk15_VVE
python python/GetDist.py chains/pk15_LSS_VVE


