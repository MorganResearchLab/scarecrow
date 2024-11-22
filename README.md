# scarecrow

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
git clone https://github.com/pachterlab/seqspec.git
scarecrow --help
scarecrow extract ./seqspec/examples/specs/SPLiT-seq/spec.yaml
```