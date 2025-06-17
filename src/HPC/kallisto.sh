#!/bin/bash
#SBATCH --job-name=kallisto
#SBATCH -o kallisto_%j.out
#SBATCH -e kallisto_%j.err
#SBATCH --ntasks=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00

# For more options see: https://slurm.schedmd.com/sbatch.html

jq=~/sharedscratch/software/jq-linux64
INDEX=
GENES=
OUT=
FASTQ=
JSON=

# Parsing named arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --index)
      INDEX="$2"
      shift
      shift
      ;;
    --genes)
      GENES="$2"
      shift
      shift
      ;;
    --out)
      OUT="$2"
      shift
      shift
      ;;
    --json)
      JSON="$2"
      shift
      shift
      ;;
    --fastq)
      FASTQ="$2"
      shift
      shift
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Compile kallisto and bustools
#module load cmake/3.23.1
#kb compile all

# Generate indexes
#kb ref -d human \
#  -i ~/sharedscratch/software/kallisto/hg38/transcriptome.idx \
#  -g ~/sharedscratch/software/kallisto/hg38/transcripts_to_genes.txt

# Run kb count to generate cell-gene count matrix
XSTR=$(${jq} -r '."kallisto-bustools"[0]."kb count" | capture("-x (?<x>[^ ]+)").x' ${JSON})
echo "kb count -t ${SLURM_CPUS_PER_TASK} -i ${INDEX} -g ${GENES} -x ${XSTR} -w NONE --h5ad --overwrite --inleave -o ${OUT} ${FASTQ}"
kb count -t ${SLURM_CPUS_PER_TASK} -i ${INDEX} -g ${GENES} -x ${XSTR} -w NONE --h5ad --overwrite --inleave -o ${OUT} ${FASTQ}

echo "Finished"
