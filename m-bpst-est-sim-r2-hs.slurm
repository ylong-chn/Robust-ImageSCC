#!/bin/sh
#SBATCH --job-name m-bpst-est-sim-r2-hs
#SBATCH --partition=normal
#SBATCH --output=/scratch/%u/%x-%A-%a.out
#SBATCH --error=/scratch/%u/%x-%A-%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ylong5@gmu.edu
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=1-00:00 

module load gnu10/10.3.0-ya
module load r/4.1.2-dx
module load python/3.10.1-5r
export PYTHONIOENCODING=utf8
 
Rscript --no-restore --quiet --no-save m-bpst-est-sim-r2-hs.R \
  --task ${SLURM_ARRAY_TASK_ID}\
  --job ${SLURM_ARRAY_JOB_ID}\