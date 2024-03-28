#!/bin/bash

# Function to generate and submit jobs
generate_and_submit_cuda() {
  MPI_NUM_TASKS=$1
  USER_CUDA=$2  # Capture the second argument as USER_CUDA
  JOB_NAME="filling_depressions$MPI_NUM_TASKS"
  OUTPUT_NAME="gflood${MPI_NUM_TASKS}"

  # Creating a temporary run script for the current number of MPI tasks
  cat > temp_run_cuda${MPI_NUM_TASKS}.sh <<EOF
#!/bin/bash
#SBATCH --time=72:30:00
#SBATCH -n ${MPI_NUM_TASKS} # MPI processes
#SBATCH -N 1 # Total number of nodes
#SBATCH --job-name=${JOB_NAME}
#SBATCH --output=${OUTPUT_NAME}.o%j
#SBATCH --gres=gpu:1 # Max of 1 GPU on Borah
#SBATCH --exclusive
#SBATCH -p gpu # Queue (partition)
module load slurm
module load mpi
mpirun -np ${MPI_NUM_TASKS} ./filling_depressions_cpucuda -F geoflood.ini --user:cuda=${USER_CUDA}
EOF

  # Submit the job
  sbatch temp_run_cuda${MPI_NUM_TASKS}.sh
}


# Loop over the desired MPI tasks
# cuda
for n in 32 16 8 4 2 1; do
  generate_and_submit_cuda $n true
done

sleep 10

generate_and_submit_cpu() {
  MPI_NUM_TASKS=$1
  USER_CUDA=$2  # Capture the second argument as USER_CUDA
  JOB_NAME="filling_depressions$MPI_NUM_TASKS"
  OUTPUT_NAME="gflood${MPI_NUM_TASKS}"

  # Creating a temporary run script for the current number of MPI tasks
  cat > temp_run_cpu${MPI_NUM_TASKS}.sh <<EOF
#!/bin/bash
#SBATCH --time=100:30:00
#SBATCH -n ${MPI_NUM_TASKS} # MPI processes
#SBATCH -N 1 # Total number of nodes
#SBATCH --job-name=${JOB_NAME}
#SBATCH --output=${OUTPUT_NAME}.o%j
# SBATCH --gres=gpu:1 # Max of 1 GPU on Borah
#SBATCH --exclusive
#SBATCH -p bsudfq # Queue (partition)
module load slurm
module load mpi
mpirun -np ${MPI_NUM_TASKS} ./filling_depressions_cpucuda -F geoflood.ini --user:cuda=${USER_CUDA}
EOF

  # Submit the job
  sbatch temp_run_cpu${MPI_NUM_TASKS}.sh
}

#cpu
for n in 32 16 8 4 2 1; do
  generate_and_submit_cpu $n false
done

# Wait for all jobs to finish (optional, remove if not needed)
#squeue -u $USER # Check if any jobs are still running
#sleep 60 # Wait for a minute (adjust based on expected job completion time)

# Function to display results
# display_results() {
#   for n in 1 2 4 8 16 32; do
#     echo "== $n processes =="
#     # Assuming %j gets replaced by job ID, and multiple files could exist
#     for file in gflood_${n}.o*; do
#       if [ -f "$file" ]; then
#         cat "$file"
#       else
#         echo "No output file found for $n processes"
#       fi
#     done
#     echo
#   done
# }

# Uncomment the following line if you wish to see the results after all jobs have finished
# display_results