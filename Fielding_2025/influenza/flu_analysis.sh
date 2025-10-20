#!/usr/bin/bash                                                                                                                                              
#SBATCH --job-name=flu_type_linked
#SBATCH --array=1-5
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=10
#SBATCH --account=hoehnlab
#SBATCH --nodelist=t01
#SBATCH --partition=preempt_t01
#SBATCH --qos=lab_priority
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=7-0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Hunter.J.Melton@dartmouth.edu
#SBATCH --exclude=t10

R CMD BATCH --no-save --no-restore flu_analysis.R ./log/flu_analysis_8_25_$SLURM_ARRAY_TASK_ID.txt
