#!/bin/bash
# Submitting batch jobs to Slurm Workload Manager

echo "Submitting run_malpasset32.qsub..."
sbatch run_malpasset32.qsub

echo "Submitting run_malpasset16.qsub..."
sbatch run_malpasset16.qsub

echo "Submitting run_malpasset8.qsub..."
sbatch run_malpasset8.qsub

echo "Submitting run_malpasset4.qsub..."
sbatch run_malpasset4.qsub

echo "Submitting run_malpasset2.qsub..."
sbatch run_malpasset2.qsub

echo "Submitting run_malpasset1.qsub..."
sbatch run_malpasset1.qsub

echo "Submitting run_malpasset64.qsub..."
sbatch run_malpasset64.qsub
