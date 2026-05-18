#!/bin/bash
# Request one node:
#SBATCH --nodes=1
#
# Specify one task:
#SBATCH --ntasks=1
#
# Number of processors for single task needed for use case (example):
#SBATCH --cpus-per-task=1
#
# Command(s) to run (example):
curr_folder="DLBCL_cell_lines"

job_submitted=(sbatch --account=${SLURM_JOB_ACCOUNT} --partition=${SLURM_JOB_PARTITION} --nodelist=${SLURM_JOB_NODELIST} --qos=${SLURM_JOB_QOS} DLBCL_cell_lines_processing.sh)

# get the job number from the command line output "Submitted batch job job_number"
job_number=$(echo "$job_submitted" | tr -dc '0-9.')

job_submitted_miseq=(sbatch --account=${SLURM_JOB_ACCOUNT} --partition=${SLURM_JOB_PARTITION} --nodelist=${SLURM_JOB_NODELIST} --qos=${SLURM_JOB_QOS} DLBCL_cell_lines_processing.sh -m)

job_miseq_number=$(echo "$job_submitted_miseq" | tr -dc '0-9.')

job_combine=(sbatch -d afterok:${job_number}:${job_miseq_number} --account=${SLURM_JOB_ACCOUNT} --partition=${SLURM_JOB_PARTITION} --nodelist=${SLURM_JOB_NODELIST} --qos=${SLURM_JOB_QOS} DLBCL_combine_miseq_neb.sh)

job_combine_number=$(echo "$job_combine" | tr -dc '0-9.')

job_define_clones=(sbatch -d afterok:${job_combine_number} --account=${SLURM_JOB_ACCOUNT} --partition=${SLURM_JOB_PARTITION} --nodelist=${SLURM_JOB_NODELIST} --qos=${SLURM_JOB_QOS} DLBCL_cell_lines_define_clones.sh)

job_define_clones_number=$(echo "$job_define_clones" | tr -dc '0-9.')

sbatch -d afterok:${job_define_clones_number} --account=${SLURM_JOB_ACCOUNT} --partition=${SLURM_JOB_PARTITION} --nodelist=${SLURM_JOB_NODELIST} --qos=${SLURM_JOB_QOS} DLBCL_cell_lines_R_analysis.sh