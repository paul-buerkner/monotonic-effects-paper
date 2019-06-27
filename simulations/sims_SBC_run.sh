#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=8G

export OMP_NUM_THREADS=1

module load r/3.5.0-python-2.7.14
export R_LIBS=$WRKDIR/R/3.5.0

srun Rscript 2018_monotonic_effects/sims_SBC_run.R
srun Rscript 2018_monotonic_effects/sims_SBC_analysis.R
