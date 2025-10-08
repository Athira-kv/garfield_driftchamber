#!/bin/bash
#SBATCH --time=10:00:00          # total expected computing time
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16       # number of CPUs
#SBATCH --output=LOG.log
#SBATCH --job-name=myjobname
#SBATCH --partition=physics
#SBATCH --mail-type=ALL
#SBATCH --mail-user=myemail

# -------- Load environment --------
module load miniconda3
source activate rootenv
# --- Run the macro ---
root -l -b -q generate_gastable.C
