#!/bin/bash
#SBATCH --job-name=samtags
#SBATCH -o samtags_%j.out
#SBATCH -e samtags_%j.err
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00

# For more options see: https://slurm.schedmd.com/sbatch.html

SAM=

# Parsing named arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --sam)
      SAM="$2"
      shift
      shift
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

samtools view --threads ${SLURM_CPUS_PER_TASK} ${SAM} | awk '{
    read_name = $1;
    tags = "";
    for (i = 12; i <= NF; i++) {
        if ($i ~ /^CB:Z:/ || $i ~ /^XM:Z:/ || $i ~ /^XP:Z:/) {
            tags = tags "\t" $i;
        }
    }
    print read_name tags;
}' > ${SAM}.tags
