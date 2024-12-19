# scarecrow

![scarecrow](img/scarecrow.png)

A toolkit for preprocessing single cell sequencing data.

### Todo

* Multi-threading is not well optimised - may be limited by I/O steps.
* Logger is currently generating two instances - could do with being a single instance.
* Currently parsing seqspec.yaml file in `scarecrow seed`, however this can probably be modified now to run without the yaml.
* Needs error handling implementing to capture missing or incorrect parameters, and unexpected file content.
* Barcode correction (`error_correct.py`) needs applying in `scarecrow harvest`.
* Summary table of barcode hits per read needs reporting in `scarecrow harvest`.

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
```bash
# Find barcode seeds
scarecrow seed <spec.yaml> <paired_fastq_R1> <paired_fastq_R2> --out <BC1_counts.csv> \
    --barcodes BC1:v1:<v1_whitelist.txt> 
scarecrow seed <spec.yaml> <paired_fastq_R1> <paired_fastq_R2> --out <BC2_counts.csv> \
    --barcodes BC2:n198:<n198_whitelist.txt>
scarecrow seed <spec.yaml> <paired_fastq_R1> <paired_fastq_R2> --out <BC3_counts.csv> \
    --barcodes BC3:v3:<v3_whitelist.txt>

# Harvest barcode positions 
scarecrow harvest <BC1_counts.csv> <BC2_counts.csv> <BC3_counts.csv> \
    --barcode_count 3 --min_distance 10 --out <barcode_positions.csv>

# Edit <barcode_positions.csv> barcode_whitelist so that each barcode (BC) has a unique name, e.g.:
#     read  start  end orientation  barcode_whitelist  read_count  read_fraction
# 0  read1     10   18     forward    [('BC1', 'v1')]        9485           0.95
# 1  read1     48   56     forward    [('BC2', 'v1')]        8393           0.84
# 2  read1     79   87     forward  [('BC3', 'n198')]        6002           0.60

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
--read2 0-100 --out cDNA.fq
```
