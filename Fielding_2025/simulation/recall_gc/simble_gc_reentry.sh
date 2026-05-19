#!/usr/bin/bash                                                                                                                                              
#SBATCH --job-name=simulate_gc_reentry
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=20
#SBATCH --account=_ACCOUNT_
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=_EMAIL_
#SBATCH --output=./log/simulate_gc_reentry.out
#SBATCH --error=./log/simulate_gc_reentry.err

module load python/3.7-Anaconda-datalad
python3 gc_reentry.py -o ./simble_sims_gc_reentry_12_18 -n 20 -p 20 --migration-rate 2