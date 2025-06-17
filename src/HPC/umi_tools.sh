#!/bin/bash
#SBATCH --job-name=umi_tools
#SBATCH -o umi_tools_%j.out
#SBATCH -e umi_tools_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00

# For more options see: https://slurm.schedmd.com/sbatch.html

module load subread/2.0.6

OUT=
BAM=
GTF=
ALIAS=

# Parsing named arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --out)
      OUT="$2"
      shift
      shift
      ;;
    --bam)
      BAM="$2"
      shift
      shift
      ;;
    --gtf)
      GTF="$2"
      shift
      shift
      ;;
    --alias)
      ALIAS="$2"
      shift
      shift
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

ID=$(basename ${BAM%.bam})

# Assigninng reads to genes
echo "Assigning reads to genes..."
echo "featureCounts -a ${GTF} -A ${ALIAS} -o ${OUT}/${ID}_gene_assigned -T ${SLURM_CPUS_PER_TASK} -R BAM ${BAM}"
featureCounts -a ${GTF} -A ${ALIAS} -o ${OUT}/${ID}_gene_assigned -T ${SLURM_CPUS_PER_TASK} -R BAM ${BAM}

echo "Sorting BAM file..."
echo "samtools sort --threads  ${SLURM_CPUS_PER_TASK} -O BAM -o ${OUT}/${ID}.featureCounts.bam ${OUT}/${ID}.bam.featureCounts.bam"
samtools sort --threads  ${SLURM_CPUS_PER_TASK} -O BAM -o ${OUT}/${ID}.featureCounts.bam ${OUT}/${ID}.bam.featureCounts.bam

echo "Indexing BAM file..."
echo "samtools index --threads  ${SLURM_CPUS_PER_TASK} ${OUT}/${ID}.featureCounts.bam"
samtools index --threads  ${SLURM_CPUS_PER_TASK} ${OUT}/${ID}.featureCounts.bam
rm ${OUT}/${ID}.bam.featureCounts.bam

#echo "Running umi_tools dedup..."
#umi_tools dedup --stdin ${OUT}/${ID}.featureCounts.bam \
#    --output-stats=${OUT}/${ID}.dedup \
#    --extract-umi-method=tag --umi-tag=UR --cell-tag=CB \
#    --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell \
#    --mapping-quality 0 --stdout ${OUT}/${ID}.featureCounts.dedup.bam

echo "Running umi_tools count..."
umi_tools count --stdin ${OUT}/${ID}.featureCounts.bam \
    --wide-format-cell-counts \
    --extract-umi-method=tag --umi-tag=UR --cell-tag=CB \
    --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell \
    --mapping-quality 0 --stdout ${OUT}/${ID}.featureCounts.counts.tsv.gz


echo "Finished"
