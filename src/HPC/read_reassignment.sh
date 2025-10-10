#!/bin/bash
#SBATCH --job-name=read_reassign
#SBATCH -o read_reassign_%j.out
#SBATCH -e read_reassign_%j.err
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00

# For more options see: https://slurm.schedmd.com/sbatch.html

TAGSA=
TAGSB=
N=100

# Parsing named arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --a)
      TAGSA="$2"
      shift
      shift
      ;;
   --b)
     TAGSB="$2"
     shift
     shift
     ;;
   --n)
     N="$2"
     shift
     shift
     ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# 1. Extract all unique barcodes from file1 if files doesn't already exist
if [ ! -f "${TAGSA}.all.barcodes" ]; then
	awk -F'\t' '{for(i=1;i<=NF;i++) if($i ~ /^CB:Z:/) {split($i,a,":"); print a[3]}}' "${TAGSA}" | sort -u > "${TAGSA}.all.barcodes"
	echo "Unique barcodes of ${TAGSA} written to ${TAGSA}.all.barcodes"
else
	echo "Unique barcodes file ${TAGSA}.all.barcodes already exists, using it"
fi

# 2. Randomly sample n barcodes
echo "Sampling ${N} barcodes from ${TAGSA}"
shuf -n "${N}" "${TAGSA}.all.barcodes" > "${TAGSA}.subset_${N}.barcodes"
if [ ! -s ${TAGSA}.subset_${N}.barcodes ]; then
	echo "Ran into a problem, ${TAGSA}.subset_${N}.barcodes is empty, exiting"
	exit 1
fi

# 3. Filter file1 on selected barcodes + extract matching read names
awk -v BARCODE_FILE="${TAGSA}.subset_${N}.barcodes" -F'\t' '
BEGIN {
    while ((getline b < BARCODE_FILE) > 0) selected[b] = 1;
}
{
    for (i = 1; i <= NF; i++) {
        if ($i ~ /^CB:Z:/) {
            split($i, arr, ":");
            bc = arr[3];
            if (bc in selected) {
                print $0;
                print $1 > "'${TAGSA}.subset_${N}_barcodes.reads'";
                break;
            }
        }
    }
}' "${TAGSA}" > "${TAGSA}.subset_${N}.tsv"

# 4. Filter file2 on read names from selected file1 lines
#grep -Ff "${TAGSA}.subset_${N}_barcodes.reads" "${TAGSB}" > "${TAGSB}.subset_${N}.tsv"
awk 'NR==FNR {a[$1]; next} $1 in a' "${TAGSA}.subset_${N}_barcodes.reads" "${TAGSB}" > "${TAGSB}.subset_${N}.tsv"

echo "Done"
echo "Filtered ${TAGSA} saved to: ${TAGSA}.subset_${N}.tsv"
echo "Filtered ${TAGSB} saved to: ${TAGSB}.subset_${N}.tsv"
