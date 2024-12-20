# scarecrow

![scarecrow](img/scarecrow.png)

A toolkit for preprocessing single cell sequencing data.

### Todo

* Multi-threading is not well optimised - may be limited by I/O steps.
* Currently parsing seqspec.yaml file in `scarecrow seed`, however this can probably be modified now to run without the yaml.
* Needs error handling implementing to capture missing or incorrect parameters, and unexpected file content.
* `scarecrow harvest` barcode counting needs implementing - need to count number of sequences with n barcodes, frequency of each barcode per barcode position, frequency of each comination of barcodes


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
scarecrow seed <spec.yaml> <paired_fastq_R1> <paired_fastq_R2> --out <BC1_counts.csv> \
    --barcodes BC1:v1:<v1_whitelist.txt> 
scarecrow seed <spec.yaml> <paired_fastq_R1> <paired_fastq_R2> --out <BC2_counts.csv> \
    --barcodes BC2:n198:<n198_whitelist.txt>
scarecrow seed <spec.yaml> <paired_fastq_R1> <paired_fastq_R2> --out <BC3_counts.csv> \
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

Finally, the `reap` tool extracts a target sequence (i.e. cDNA) and its associated quality values from either `--read1` or `--read2`. This sequence data is written to a new fastq file, and the combination of cell barcodes identified near the positions previously estimated are appened to the sequence header. The `jitter` value is the flanking distance to extend the barcode search for from the start and end positions of the barcode peaks. `mismatch` is the maximum number of mismatches between an expected and observed barcode sequence. *Barcode error correction remains to be implemented*.

```bash
# Reap target sequence from fastqs (TBC)
scarecrow reap <paired_fastq_R1> <paired_fastq_R2> --barcode_positions <barcode_positions.csv> \
    --barcodes BC1:v1:<v1_whitelist.txt> BC2:v1:<v1_whitelist.txt> BC3:n198:<n198_whitelist.txt> \
    --jitter 5 --mismatches 1 --read2 0-100 --out cDNA.fastq
```

### Example
```bash
# seed
scarecrow seed /Users/s14dw4/Documents/Repos/scarecrow/specs/evercode/evercode-v3.yaml \
    /Users/s14dw4/Documents/Repos/scarecrow/specs/evercode/R1.fastq \
    /Users/s14dw4/Documents/Repos/scarecrow/specs/evercode/R2.fastq \
    --barcodes BC1:v1:/Users/s14dw4/Documents/Repos/scarecrow/specs/evercode/BC1.txt \
    --out ./results/barcodes_BC1.csv

# harvest
FILES=(./results/barcodes*csv)
scarecrow harvest ${FILES[@]} --barcode_count 3 --min_distance 10 \
    --out ./barcode_positions.csv

# reap
scarecrow reap /Users/s14dw4/Documents/Repos/scarecrow/specs/evercode/R1.fastq \
    /Users/s14dw4/Documents/Repos/scarecrow/specs/evercode/R2.fastq \
    -p ./barcode_positions.csv -j 5 -m 1 \
    --barcodes \
        BC1:v1:/Users/s14dw4/Documents/Repos/scarecrow/barcodes/bc_data_v1.barcodes \
        BC2:v1:/Users/s14dw4/Documents/Repos/scarecrow/barcodes/bc_data_v1.barcodes \
        BC3:n198:/Users/s14dw4/Documents/Repos/scarecrow/barcodes/bc_data_n198_v5.barcodes \
    --read2 0-100 --out ./cDNA.fq
```
