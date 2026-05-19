#!/usr/bin/bash 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=6gb
#SBATCH --time=21-00:00
#SBATCH --array=1-4

final_folder=output_5_2/sel/
mkdir -p $final_folder

# exit if this line fails
Rscript tyche_reentry_test_constraints.R sel 5_2 || exit 1

cp -r /scratch/5_2_sel/* $final_folder