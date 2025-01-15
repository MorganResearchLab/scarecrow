# scarecrow

![scarecrow](img/scarecrow.png)

A toolkit for preprocessing single cell sequencing data.

### Todo

* Needs error handling implementing to capture missing or incorrect parameters, and unexpected file content.
* `scarecrow reap` needs checking after optimizing for efficiency
* Replace --read1 --read2 flags in reap with --extract <read>:<range>
* Consider adding --udi <read>:<range> to reap
* Refactor harvest for efficiency
* Increase font size on harvest plot


## Environment
```bash
mamba create --name scarecrow python=3.12
mamba activate scarecrow
mamba install git pip
pip install git+https://github.com/pachterlab/seqspec.git
pip install git+https://github.com/MorganResearchLab/scarecrow.git
mamba install numpy, scipy, pandas, seaborn
```

## Workflow
There are currently three tools: `seed`, `harvest`, `reap`, each of which feeds into the next. Using a small subset of paired-end sequencing reads (~10K ), the `seed` tool identifies barcode seeds within reads and records the counts to a CSV file. This should be run on each barcode whitelist if the library structure is uncertain.

```bash
# Find barcode seeds
scarecrow seed --fastqs <paired_fastq_R1> <paired_fastq_R2> --strands pos neg --out <BC1_counts.csv> \
    --barcodes BC1:v1:<v1_whitelist.txt> 
scarecrow seed --fastqs <paired_fastq_R1> <paired_fastq_R2> --strands pos neg --out <BC2_counts.csv> \
    --barcodes BC2:n198:<n198_whitelist.txt>
scarecrow seed --fastqs <paired_fastq_R1> <paired_fastq_R2> --strands pos neg --out <BC3_counts.csv> \
    --barcodes BC3:v3:<v3_whitelist.txt>
```

The resulting CSV files are then processed with `harvest` to identify barcode alignment peaks. The minimum distance between barcodes can be provided, and the number of peaks (--barcode_count) to return. The results are recorded to a CSV file and a histogram of the peaks is plotted to a PNG file. The CSV table should be edited to ensure that the barcode_whitelist values are unique.

```bash
# Harvest barcode positions 
scarecrow harvest <BC1_counts.csv> <BC2_counts.csv> <BC3_counts.csv> \
    --barcode_count 3 --min_distance 10 --out <barcode_positions.csv>

# Edit <barcode_positions.csv> barcode_whitelist so that each barcode (BC) has a unique name, e.g.:
#     read  start  end orientation  barcode_whitelist  read_count  read_fraction
# 0  read1     10   18     forward    [('BC1', 'v1')]        9485           0.95
# 1  read1     48   56     forward    [('BC2', 'v1')]        8393           0.84
# 2  read1     79   87     forward  [('BC3', 'n198')]        6002           0.60
```

Finally, the `reap` tool extracts a target sequence (i.e. cDNA) and its associated quality values from either `--read1` or `--read2`. This sequence data is written to a new fastq file, and the combination of cell barcodes identified near the positions previously estimated are appened to the sequence header. The `jitter` value is the flanking distance to extend the barcode search for from the start and end positions of the barcode peaks. `mismatch` is the maximum number of mismatches between an expected and observed barcode sequence.

```bash
# Reap target sequence from fastqs (TBC)
scarecrow reap --fastqs <paired_fastq_R1> <paired_fastq_R2> --barcode_positions <barcode_positions.csv> \
    --barcodes BC1:v1:<v1_whitelist.txt> BC2:v1:<v1_whitelist.txt> BC3:n198:<n198_whitelist.txt> \
    --jitter 5 --mismatches 1 --read2 0-100 --out cDNA.fastq.gz
```

### Parse Evercode WTv3 Example
```bash
OUTDIR=/uoa/home/s14dw4/scarecrow_test/MorganLab

# seed
R1=/uoa/scratch/shared/Morgan_Lab/Parse/LV6000778329-SO8229-University-of-Aberdeen-Morgan-Lab-Lib1-2024-12-10_S1_L001_R1_001.fastq.gz
R2=/uoa/scratch/shared/Morgan_Lab/Parse/LV6000778329-SO8229-University-of-Aberdeen-Morgan-Lab-Lib1-2024-12-10_S1_L001_R2_001.fastq.gz

# barcodes
BARCODES=(
  "BC1:R1_v3:${OUTDIR}/barcodes/bc_data_n123_R1_v3_5.csv"
  "BC2:v1:${OUTDIR}/barcodes/bc_data_v1.csv"
  "BC3:R3_v3:${OUTDIR}/barcodes/bc_data_R3_v3.csv"
  )

# Exract 100K reads
zcat -c ${R1} | head --l 400000 > ${OUTDIR}/fastq/100K_1.fastq 
zcat -c ${R2} | head --l 400000 > ${OUTDIR}/fastq/100K_2.fastq 

# Identify barcode seeds
for BARCODE in ${BARCODES[@]}
do
    scarecrow seed --fastqs ${R1} ${R2} --strands pos neg \
      -o ${OUTDIR}/results/barcodes_${BARCODE%:*:*}.csv --barcodes ${BARCODE}
done    

# Harvest barcode peaks
FILES=(${OUTDIR}/results/barcodes_*csv)
scarecrow harvest ${FILES[@]} --barcode_count 3 --min_distance 11 \
    --out ${OUTDIR}/barcode_positions.csv

# Reap cDNA
scarecrow reap --fastqs ${R1} ${R2} -p ${OUTDIR}/barcode_positions.csv \
    -j 5 -m 1 --barcodes ${BARCODES[@]} --read1 0-64 --out ./cDNA.fq.gz  
```

To process 150M read pairs with `scarecrow reap` takes around 90 mins, requires 1 core and ~ 200 MB RAM.



# Testing on laptop
```bash
R1=100K_1.fastq
R2=100K_2.fastq
BARCODES=("BC1:R1_v3:/Users/s14dw4/Documents/scarecrow_test/barcodes/bc_data_n123_R1_v3_5.barcodes"
    "BC2:v1:/Users/s14dw4/Documents/scarecrow_test/barcodes/bc_data_v1.barcodes"
    "BC3:R3_v3:/Users/s14dw4/Documents/scarecrow_test/barcodes/bc_data_R3_v3.barcodes")
for BARCODE in ${BARCODES[@]}
do
    scarecrow seed --fastqs ${R1} ${R2} --strands pos neg \
      -o ./results/barcodes_${BARCODE%:*:*}.csv --barcodes ${BARCODE}
done

FILES=(./results/barcodes_BC*csv)
scarecrow harvest ${FILES[@]} --barcode_count 3 --min_distance 11 --out barcode_positions.csv
```
