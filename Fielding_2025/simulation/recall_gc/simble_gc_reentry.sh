#!/usr/bin/bash                                                                                                                                              
#SBATCH --job-name=simulate_gc_reentry
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=10
#SBATCH --account=hoehnlab-share
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Hunter.J.Melton@dartmouth.edu
#SBATCH --exclude=t10
#SBATCH --output=./log/simulate_gc_reentry.out
#SBATCH --error=./log/simulate_gc_reentry.err

module load python/3.7-Anaconda-datalad
conda activate DiscoveryPythonEnv
python3 gc_reentry.py -o ./simble_sims_gc_reentry_8_28 -n 20 -p 10 --migration-rate 2