#!/usr/bin/bash                                                                                                                                              
#SBATCH --job-name=flu_4_13
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --account=_ACCOUNT_
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=14-0
#SBATCH --array=1-3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=_EMAIL_


R CMD BATCH --no-save --no-restore flu_analysis.R ./log/flu_analysis_4_13_$SLURM_ARRAY_TASK_ID.txt
