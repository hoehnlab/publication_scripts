#!/usr/bin/bash                                                                                                                                              
#SBATCH --job-name=gc_neu_reentry_sims
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --account=hoehnlab-share
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=21-0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Hunter.J.Melton@dartmouth.edu
#SBATCH --exclude=t10

R CMD BATCH --no-save --no-restore gc_reentry_sims_uniform_neutral.R ./log_un/gc_reentry_sims_uniform_neutral_9_9_$SLURM_ARRAY_TASK_ID.txt
