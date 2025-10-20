#!/bin/bash
# Job name:
#SBATCH --job-name=simble-validation-final
#
# Request one node:
#SBATCH --array=1-7
#
# Specify one task:
#SBATCH --ntasks=1
#
# Number of processors for single task needed for use case (example):
#SBATCH --cpus-per-task=20
#
# Wall clock limit:
#SBATCH --time=10:00:00
#
## Command(s) to run (example):
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate sim


case $SLURM_ARRAY_TASK_ID in
1)
    simble -o ~/lab/jessie/simble-validation/200gen_neutral866228 --seed 866228 --neutral -n 100 -p 20 -s 0 200 50 --sample-size 20 -q
    ;;
2)
    simble -o ~/lab/jessie/simble-validation/200gen_uniform358339 --seed 358339 --neutral --uniform -n 100 -p 20 -s 0 200 50 --sample-size 20 -q
    ;;
3)
    simble -o ~/lab/jessie/simble-validation/200gen_selection233853 --seed 233853 -n 100 -p 20 -s 0 200 50 --sample-size 20 -q
    ;;
4)
    simble -o ~/lab/jessie/simble-validation/150gen_selection_100clones_differentiation_study_seed344276 --seed 344276 -n 100 -p 20 -s 0 150 150 --sample-size 0 --sample-size-other 1000 --migration 5 -q
    ;;
5)
    python ~/simble_addons/cross_reactivity_test.py --seed 226348 -o ~/lab/jessie/simble-validation/cross_reactivity_226348/ -n 10 -p 10 -s 0 150 150 --sample-size 0 --sample-size-other 1000 --migration-rate 5 -q
    ;;
6)
    simble -o ~/lab/jessie/simble-validation/500gen_selection_150clones_seed253437 --seed 253437 -n 150 -p 20 -s 0 500 50 --sample-size 50 -q
    ;;
7)
    python ~/simble_addons/mutations_per_site.py --seed 754659 -o ~/lab/jessie/simble-validation/mutations_per_site_754659/ -s 0 200 50 --sample-size 20 -q

esac