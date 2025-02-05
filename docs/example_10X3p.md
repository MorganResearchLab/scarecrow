<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[Back to root](root.md)

# Example: 10X 3'

<img style="float:right;width:100%;" src="../img/10X3p.svg" alt="scarecrow"/>

Stay organised - create a folder for the project to keep things tidy.

```bash
PROJECT=./scarecrow/examples/10X3p
mkdir -p ${PROJECT}
```

Download 10X 3' data from [https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1106903](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1106903). Instead of `wget`, this time we use `fasterq-dump` as we need one of the technical reads. These accessions are data for two technical replicates.

```bash
mkdir -p ${PROJECT}/fastq
ACCS=(SRR28867562 SRR28867563)
for ACC in ${ACCS[@]}
do
    prefetch --output-directory ${PROJECT}/fastq ${ACC}
    fasterq-dump ${PROJECT}/fastq/${ACC} -e 2 --split-files --include-technical --force --outdir ${PROJECT}/fastq
    gzip ${PROJECT}/fastq/${ACC}_1.fastq
    gzip ${PROJECT}/fastq/${ACC}_2.fastq
    gzip ${PROJECT}/fastq/${ACC}_3.fastq
    gzip ${PROJECT}/fastq/${ACC}_4.fastq
done
```

### 1. Extract subset of 1M reads for profiling

```bash
mkdir -p ${PROJECT}/fastq/subset
FILES=(${PROJECT}/fastq/*.fastq.gz)
for FILE in ${FILES[@]}
do
    ID=$(basename ${FILE%.fast*})
    zcat -c ${FILE} | head --l 4000000 | gzip > ${PROJECT}/fastq/subset/${ID}.fastq.gz
done
```

### 2. Generate barcode match profiles

Chromium barcode whitelists for different chemistry versions are available at [https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html). The library for this sample used the v3.1 chemistry, so the barcode whitelist to use is 3M-february-2018.txt.gz. The file should be downloaded to `${PROJECT}/barcodes`, unzipped and checked that it contains a single barcode per line with no headers. 

We can now run `scarecrow seed` to process the barcode whitelist. The below example is for a SLURM HPC, but will work on a standard PC by omitting the `sbatch` line.

***There are significantly more barcodes with 10X than Parse. May need to refactor seed to work with a set like reap. Also check if subset of barcodes returns equivalent profile.***

```bash
mkdir -p ${PROJECT}/barcode_profiles
THREADS=1
FASTQS=(${PROJECT}/fastq/subset/*.fastq.gz)
BARCODE=BC1:3M-FEB-2018:${PROJECT}/barcodes/3M-february-2018.txt
sbatch --ntasks ${THREADS} --mem 48G --time=01:00:00 -o seed.%j.out -e seed.%j.err \
    scarecrow seed \
        --threads ${THREADS} \
        --fastqs ${FASTQS[@]} --strands pos neg \
        --barcodes ${BARCODE} \
        --out ${PROJECT}/barcode_profiles/barcodes.${BARCODE%%:*}.csv
```