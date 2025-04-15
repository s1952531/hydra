#!/bin/csh

#SBATCH --ntasks-per-node=1
#SBATCH --account=ec12
#SBATCH --job-name=hydra-basin
#SBATCH --output=log
#SBATCH --mem-per-cpu=4G
#SBATCH --time=5-00:00:00

caps
