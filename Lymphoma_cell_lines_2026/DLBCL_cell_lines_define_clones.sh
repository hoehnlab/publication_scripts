#!/bin/bash
# Request one node:
#SBATCH --array=1-16
#
# Specify one task:
#SBATCH --ntasks=1
#
# Number of processors for single task needed for use case (example):
#SBATCH --cpus-per-task=8
#
# Wall clock limit:
#SBATCH --time=05:00:00
#
## Command(s) to run (example):
export IGDATA=~/test-igblast
curr_folder=combined_igblast


dodefineclones() {
    local name=$1
    DefineClones.py -d ${curr_folder}/${name}.tsv --format airr --act set --model ham --norm len --sym min --dist 0.1 --nproc 8 --outdir define_clones_combined
}

# List of cell lines to process
cell_lines=(
    "Bcl2"
    "CBP1"
    "CBP2"
    "CMBP3"
    "CMBP7"
    "EBMP2"
    "EBMP3"
    "EZBP1"
    "EZBP6"
    "MCDP2"
    "MCDP4"
    "MCDP6"
    "SN2F"
    "SN2M"
    "WP2"
    "WP6"
)

n=${SLURM_ARRAY_TASK_ID}-1


cell_line=${cell_lines[$n]}


echo "Processing cell line: $cell_line"

dodefineclones ${cell_line}-final_collapse-unique_atleast-2_igblast
