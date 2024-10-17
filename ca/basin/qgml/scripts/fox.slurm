#!/bin/bash

##SBATCH --exclusive
#SBATCH --ntasks-per-node=1
#SBATCH --account=ec12
#SBATCH --job-name=qgbasin
#SBATCH --output=log
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0-01:00:00
##SBATCH --partition=bigmem

# Recommended safety settings
set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

rm log
./caps
