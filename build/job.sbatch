#!/bin/bash
#SBATCH --time=02:00:10          # total expected computing time
#SBATCH --ntasks=1
#SBATCH --nodes=1                # number of nodes
#SBATCH --cpus-per-task=16       # number of CPUs
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=16G
#SBATCH --output=LOG.log
#SBATCH --job-name=myjobname
#SBATCH --partition=secondary    # <-- use this one
#SBATCH --mail-type=ALL
#SBATCH --mail-user=myemail

# -------- Load environment --------
module load miniconda3
source activate rootenv

# -------- Run your executable --------
./Gem
