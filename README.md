# scarecrow

![scarecrow](img/scarecrow.png)

A toolkit for preprocessing single cell sequencing data.

### Issues

* Multi-threading is not well optimised - may be limited by I/O steps.
* Logger is currently generating two instances - could do with being a single instance.
* Currently parsing seqspec.yaml file in `scarecrow seed`, however this can probably be modified now to run without the yaml.

## Environment
```bash
mamba create --name scarecrow python=3.12
mamba activate scarecrow
mamba install git pip
pip install git+https://github.com/pachterlab/seqspec.git
pip install git+https://github.com/MorganResearchLab/scarecrow.git
mamba install numpy, scipy, pandas
```

## Workflow
```bash
# Find barcode seeds
scarecrow seed <spec.yaml> <fastqs> --barcodes BC1:<barcode_whitelist.txt> --out <BC1_counts.csv>

# Harvest barcode positions 
scarecrow harvest <BCx_counts> --barcode_count 3 --min_distance 10

# Reap target sequence from fastqs (TBC)
scarecrow reap
```

### Example
```bash
# seed
scarecrow seed ./scarecrow/specs/evercode/evercode-v3.yaml \
/Users/s14dw4/Documents/Repos/scarecrow/specs/evercode/R1.fastq \
/Users/s14dw4/Documents/Repos/scarecrow/specs/evercode/R2.fastq \
--barcodes BC1:/Users/s14dw4/Documents/Repos/scarecrow/specs/evercode/BC1.txt

# harvest
scarecrow harvest ./barcode_counts.csv --barcode_count 3 --min_distance 10
```