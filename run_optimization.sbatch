#!/bin/bash
#SBATCH --qos=open
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH --time=47:59:00
#SBATCH --job-name=csdnoise-opt
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

cd $SLURM_SUBMIT_DIR

echo "Job started on $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Threads: $SLURM_CPUS_PER_TASK"

JULIA_VERSION=$(julia --version | grep -oP '\d+\.\d+\.\d+')
REQUIRED_VERSION=$(grep '^julia_version' Manifest.toml | grep -oP '\d+\.\d+\.\d+')

echo "Julia version: $JULIA_VERSION"
echo "Required version: $REQUIRED_VERSION"

if [ "$JULIA_VERSION" != "$REQUIRED_VERSION" ]; then
    echo "ERROR: Julia version mismatch!"
    echo "Expected: $REQUIRED_VERSION"
    echo "Found: $JULIA_VERSION"
    echo "Cancelling job..."
    scancel $SLURM_JOB_ID
    exit 1
fi

echo "Julia version check passed. Starting optimization..."

julia --threads=$SLURM_CPUS_PER_TASK --startup-file=no --history-file=no --project=CSDNoiseCore ./scripts/gridsearch-optimization.jl

echo "Job completed on $(date)"

sacct -j $SLURM_JOB_ID --format=JobID,JobName,MaxRSS,Elapsed,TotalCPU,State
