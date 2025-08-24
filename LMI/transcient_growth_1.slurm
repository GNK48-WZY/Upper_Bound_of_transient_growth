#!/bin/bash
### partition 1: RM [default,use one or more full nodes],
### partition 2: RM-shared [use only part of one node] 128 cpu per node
# partition 3: RM-512 [use one or more full 512GB-memory nodes]
###SBATCH -p RM
#SBATCH -p EM
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=24:00:00
###SBATCH --mail-user=zhengyang.wei@uconn.edu      # Destination email address
###SBATCH --mail-type=ALL                       # Event(s) that triggers email notification
#SBATCH --job-name=matlab_job
#SBATCH --output=matlab_output_%j
export SLURM_EXPORT_ENV=ALL
#export I_MPI_FABRICS=shm,tcp

##the slurm number to restart simulation... This need full state to be stored.
SUBMITDIR=$SLURM_SUBMIT_DIR
WORKDIR=/ocean/projects/phy240052p/zwei4/testing/matlab_$SLURM_JOB_ID
# mkdir -p "$WORKDIR" && cp -r *.m "$WORKDIR" && cp submit_dedalus_bridges2 "$WORKDIR" && cd "$WORKDIR" || exit -1
mkdir -p "$WORKDIR" && cp -r *.m "$WORKDIR" && cd "$WORKDIR" || exit -1

# submit with base envs from local computer
###source activate base
###conda activate dedalus3
###mpiexec -n $SLURM_NTASKS python3 RBC.py

###cd "$SUBMITDIR" && cp dedalus_output_$SLURM_JOB_ID "$WORKDIR"
module load matlab/R2022b
matlab -nodisplay -nosplash -nodesktop -r "transcient_growth"\;exit\;
#matlab -nodisplay -nosplash -nodesktop -logfile transcient_growth.log < transcient_growth.m > transcient_growth.out

