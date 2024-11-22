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
```bash
git clone https://github.com/cellatlas/cellatlas.git # Cloned for SPLiTSeq example fastq files
git clone https://github.com/pachterlab/seqspec.git  # Cloned for SPLiTSeq v0.3.0 spec.yaml
scarecrow --help
scarecrow extract ./seqspec/examples/specs/SPLiT-seq/spec.yaml \
    ./cellatlas/examples/rna-splitseq/fastqs/R1.fastq.gz \
    ./cellatlas/examples/rna-splitseq/fastqs/R2.fastq.gz \
    -vp
```
