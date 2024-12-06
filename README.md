# scarecrow

![scarecrow](img/scarecrow.png)

A toolkit to parse [seqspec](https://github.com/pachterlab/seqspec) `yaml` files for processing of single cell sequencing data.

### Issues

* Multi-threading is not well optimised - may be limited by I/O steps
* Logger is currently generating two instances - could do with being a single instance

## Environment
```bash
mamba create --name scarecrow python=3.12
mamba activate scarecrow
mamba install git pip
pip install git+https://github.com/pachterlab/seqspec.git
pip install git+https://github.com/MorganResearchLab/scarecrow.git
```

#### Test: Evercode
Currently, the `read_id` in `spec.yaml` must match the full path to the fastq files, and fastq files must also be passed in full to `scarecrow`. Within the specs folder for the repo are some example cases. For Parse Evercode WTv2 there is a draft `spec.yaml` (*element positions currently incorrect*), along with barcode details for BC1, BC2 and BC3 (*there is repitition of barcodes across the two whitelists). Also included are the i5 and i7 whitelists for reference. The two fastq files are the first 100 reads from `ERR12167395` (*original fastq is 11.4M read pairs*).

`scarecrow barcodes` will search each read pair for the barcodes on each whitelist and return a CSV file of hits (read, read pair number, barcode, start, end, mismatches) within 1 mismatch each barcode.

`scarecrow extract` will extract a target sequence element (-t) annotated in `spec.yaml`, and record the barcode elements (-r) in the sequence header, returning a single fastq file (header, sequence + qualities).

```bash
# Search for barcodes (also works for any sequence motif on a whitelist, e.g. i5, i7)
scarecrow barcodes ./scarecrow/specs/evercode/evercode-v3.yaml \
    /Users/s14dw4/Documents/Repos/scarecrow/specs/evercode/R1.fastq \
    /Users/s14dw4/Documents/Repos/scarecrow/specs/evercode/R2.fastq \
    --barcodes \
    BC1:/Users/s14dw4/Documents/Repos/scarecrow/specs/evercode/BC1.txt \
    BC2:/Users/s14dw4/Documents/Repos/scarecrow/specs/evercode/BC2-3.txt \
    BC3:/Users/s14dw4/Documents/Repos/scarecrow/specs/evercode/BC2-3.txt

# Extract cDNA element
scarecrow extract ./scarecrow/specs/evercode/evercode-v3.yaml \
    /Users/s14dw4/Documents/Repos/scarecrow/specs/evercode/R1.fastq \
    /Users/s14dw4/Documents/Repos/scarecrow/specs/evercode/R2.fastq \
    -o ./cDNA.fq -t cdna -r BC1 BC2 BC3
```

