#!/bin/bash
#SBATCH --job-name= samplig_Recon3D_ACHR
#SBATCH --mail-user=deepen.data@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err
#SBATCH --partition=slims
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -L matlab,matlab-distrib_computing
#SBATCH --mem=48000

module load Matlab/2017
module load gurobi/7.5.1

mkdir -p /tmp/${SLURM_JOB_ID}/{lcj,mcr}
export MCR_CACHE_ROOT=/tmp/${SLURM_JOB_ID}/mcr

export MATLAB_N=$1

matlab -nodisplay -nosplash -nodesktop < /home/aacevedo/NLHPC/Scripts/Leftraru_shrinking_Recon2_ACHR.m
