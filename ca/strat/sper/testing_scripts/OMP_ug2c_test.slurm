#!/bin/bash
#SBATCH --job-name=OMP_test_ug2c_straka # Name of the job
#SBATCH --output=OMP_ug2c_test_%j.log   # Stdout will go to log
#SBATCH --error=OMP_ug2c_test_%j.err    # Stderr will go to caps_<jobid>.err
#SBATCH --ntasks=1                      # Number of MPI tasks (processes)
#SBATCH --cpus-per-task=48              # Number of CPU cores for this task
#SBATCH --time=03:00:00                 # Maximum runtime (HH:MM:SS)
#SBATCH --partition=normal              # Partition to submit to
#SBATCH --nvram-options=1LM:1000

#   We want one thread per physical core
export OMP_PLACES=cores
#   Bind threads to cores in order
export OMP_PROC_BIND=close

#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=2

export OMP_STACKSIZE=1G
ulimit -s unlimited

echo "=== OpenMP Settings ==="
echo "OMP_NUM_THREADS = $OMP_NUM_THREADS"
echo "OMP_STACKSIZE  = $OMP_STACKSIZE bytes"
echo "OMP_PLACES    = $OMP_PLACES"
echo "OMP_PROC_BIND = $OMP_PROC_BIND"
echo "======================="

# Run the application
srun --hint=nomultithread /home/nx01/nx01/a_holmes/hydra/ca/strat/sper/testing_scripts/ugrid2con_test.sh
