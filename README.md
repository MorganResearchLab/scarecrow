# scarecrow

![scarecrow](img/scarecrow.png)

A toolkit to parse [seqspec](https://github.com/pachterlab/seqspec) `yaml` files for processing of single cell sequencing data.

## Environment
```bash
mamba create --name scarecrow python=3.12
mamba activate scarecrow
mamba install git pip
pip install git+https://github.com/pachterlab/seqspec.git
pip install git+https://github.com/MorganResearchLab/scarecrow.git
```

## Test
The fastq files must match the `read_id` in the spec.yaml. A modified file is provided with scarecrow (`./scarecrow/specs/splitseq/splitseq.yaml`), however the `read_id` requires updating on a user-by-user basis if the fastq files are not stored in a shared space.

```bash
git clone https://github.com/cellatlas/cellatlas.git # Cloned for SPLiTSeq example fastq files
git clone https://github.com/pachterlab/seqspec.git  # Cloned for SPLiTSeq v0.3.0 spec.yaml
scarecrow --help
scarecrow extract ./scarecrow/specs/splitseq/splitseq.yaml \
    /Users/s14dw4/Documents/Repos/cellatlas/examples/rna-splitseq/fastqs/R1.fastq.gz \
    /Users/s14dw4/Documents/Repos/cellatlas/examples/rna-splitseq/fastqs/R2.fastq.gz \
    -o ./cDNA.fq -r UMI Round_1_BC Round_2_BC Round_3_BC 
```

Expected console output:

```bash
seqspec print ./scarecrow/specs/splitseq/splitseq.yaml

                                          ┌─'P5:29'
                                          ├─'Spacer:8'
                                          ├─'Read_1_primer:33'
                                          ├─'cDNA:100'
                                          ├─'RT_primer:15'
                                          ├─'Round_1_BC:8'
                                          ├─'linker_1:30'
──────────────────── ──rna────────────────┤
                                          ├─'Round_2_BC:8'
                                          ├─'Linker_2:30'
                                          ├─'Round_3_BC:8'
                                          ├─'UMI:10'
                                          ├─'Read_2_primer:22'
                                          ├─'Round_4_BC:6'
                                          └─'P7:24'

Library elements identified by seqspec.get_index_by_primer

/Users/s14dw4/Documents/Repos/cellatlas/examples/rna-splitseq/fastqs/R1.fastq.gz
        cDNA: 0-100
        RT_primer: 100-115
        Round_1_BC: 115-123
        linker_1: 123-140

/Users/s14dw4/Documents/Repos/cellatlas/examples/rna-splitseq/fastqs/R2.fastq.gz
        UMI: 0-10
        Round_3_BC: 10-18
        Linker_2: 18-48
        Round_2_BC: 48-56
        linker_1: 56-86

Processing cDNA and writing to ./cDNA.fq
Processing Batches: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████| 1000/1000 [00:05<00:00, 187.04it/s]

Total read pairs processed: 1000000

Barcode Distribution
Bin Ranges (Read Counts) | Histogram
1-2                  | #################### (740165)
2-4                  | ##                   (82638)
4-8                  |                      (6911)
8-16                 |                      (2030)
16-31                |                      (423)
31-62                |                      (60)
62-124               |                      (40)
124-247              |                      (47)
247-491              |                      (15)
491-977              |                      (1)

Total unique barcodes: 832330
Min count: 1
Max count: 977

kallisto bus -x 0,115,123,1,10,18,1,48,56:1,0,10:0,0,100
```

Output files should include `cDNA.rq`, `cDNA.fq.barcode_counts.csv`, and `scarecrow.log`. The barcode counts in `cDNA.fq.barcode_counts.csv` show incorrect barcodes (top 5 listed below):

```bash
umi_barcodes,Count
AGATTCGTCA_ACGATCAG_GCTACGCT_CCATAGTT,977
CATTCCTAGG_AAAAAAAA_GCTACGCT_TGTCAGAG,461
CGTTACTAGG_AAAAAAAA_GCTACGCT_TGTCAGAG,434
TTTAACCCGG_AAAAAAAA_GCTACGCT_TGTCAGAG,387
CGTACTCTGG_AAAAAAAA_GCTACGCT_TGTCAGAG,360
```

The barcode whitelists for the `seqspec` SPLiT-seq example are available in the repo (`./seqspec/examples/specs/SPLiT-seq/`).


## Check barcodes
To check if barcode sequences from a whitelist are present in a fastq file.

```bash
python ./scarecrow/src/scarecrow/barcode_check.py \
    ./seqspec/examples/specs/SPLiT-seq/onlist_round1.txt \
    ./cellatlas/examples/rna-splitseq/fastqs/R1.fastq.gz \
    ./Round_1_BC_counts.txt
```

