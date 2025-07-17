#!/bin/sh
#SBATCH --job-name m-bpst-scc-sim
#SBATCH --partition=normal
#SBATCH --output=/scratch/%u/%x-%A-%a.out 
#SBATCH --error=/scratch/%u/%x-%A-%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4GB 
#SBATCH --time=4-00:00 

## Load the relevant modules needed for the job
module load gnu10/10.3.0-ya
module load r/4.1.2-dx
 
## Start the job
Rscript --no-restore --quiet --no-save m-bpst-scc-sim.R \
  --task ${SLURM_ARRAY_TASK_ID}\
  --job ${SLURM_ARRAY_JOB_ID}\
