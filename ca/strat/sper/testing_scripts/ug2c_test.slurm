#!/bin/bash
#SBATCH --job-name=test_ug2c_straka    # Name of the job
#SBATCH --output=ug2c_test_%j.log                    # Stdout will go to log
#SBATCH --error=ug2c_test_%j.err             # Stderr will go to caps_<jobid>.err
#SBATCH --ntasks=1                      # Number of MPI tasks (processes)
#SBATCH --cpus-per-task=48              # Number of CPU cores for this task
#SBATCH --time=03:00:00                 # Maximum runtime (HH:MM:SS)
#SBATCH --partition=normal              # Partition to submit to
#SBATCH --nvram-options=1LM:1000

# Run the application
srun /home/nx01/nx01/a_holmes/hydra/ca/strat/sper/testing_scripts/ugrid2con_test.sh
