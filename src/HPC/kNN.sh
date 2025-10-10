#!/bin/bash
#SBATCH --job-name=kNN
#SBATCH -o kNN_%j.out
#SBATCH -e kNN_%j.err
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00

# For more options see: https://slurm.schedmd.com/sbatch.html

kNN=./scarecrow/scripts/kNN.R
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

echo "${kNN} ${IN}"
${kNN} ${IN}
