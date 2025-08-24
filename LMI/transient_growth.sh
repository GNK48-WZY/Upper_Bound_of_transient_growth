#!/bin/bash -l
#SBATCH --job-name=matlab
#SBATCH --time=12:00:00
###SBATCH --partition=lo-core
###SBATCH --partition=lo-core # This can be as long as 7 days
###SBATCH --account=chl23026
#SBATCH --partition=priority
###SBATCH --qos=me_epyc
###SBATCH --partition=priority # This can run infinite time with priority in queue
#SBATCH --mem=200G
###SBATCH --no-requeue
###SBATCH --mem-per-cpu=4090M
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=matlab_history_%j
#---------------------------------------------------------------------
# SLURM job script to run serial MATLAB
#---------------------------------------------------------------------

SUBMITDIR=$SLURM_SUBMIT_DIR
WORKDIR=/scratch/chl23026/zhw22003/matlab_$SLURM_JOB_ID
mkdir -p "$WORKDIR" && cp -r *.m "$WORKDIR" && cd "$WORKDIR" || exit -1

module load matlab/R2023b
matlab -nodisplay -nosplash -nodesktop -r "transcient_growth_new"\;exit\;

