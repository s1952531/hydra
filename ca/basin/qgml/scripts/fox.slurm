#!/bin/csh

#SBATCH --ntasks-per-node=1
#SBATCH --account=ec12
#SBATCH --job-name=-hydra
#SBATCH --output=log
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0-18:00:00

caps
