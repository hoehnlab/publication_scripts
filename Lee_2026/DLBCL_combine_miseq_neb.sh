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
# Wall clock limit:
#SBATCH --time=02:00:00
#
# Command(s) to run (example):

# run the R script for analysis
Rscript DLBCL_combine_miseq_neb.R