#!/bin/bash
#SBATCH --job-name=counts2mtx
#SBATCH -o counts2mtx_%j.out
#SBATCH -e counts2mtx_%j.err
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00

# For more options see: https://slurm.schedmd.com/sbatch.html

counts2mtx=./scarecrow/src/HPC/counts2mtx.R
IN=


# Parsing named arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --in)
      IN="$2"
      shift
      shift
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

echo "${counts2mtx} ${IN}"
${counts2mtx} ${IN}
