#!/usr/bin/bash 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=6gb
#SBATCH --time=21-00:00
#SBATCH --array=1-4

final_folder=output_4_13/un/
mkdir -p $final_folder

# exit if this line fails
Rscript tyche_gc_reentry.R un 4_13 || exit 1

cp -r /scratch/4_13_un/* $final_folder