<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[Back to root](root.md)

# Example: Parse Evercode WTv2

Stay organised - create a folder for the project to keep things tidy.

```bash
PROJECT=./scarecrow/examples/evercode
mkdir -p ${PROJECT}
```

Download Evercode WTv2 data from [https://www.ebi.ac.uk/ena/browser/view/PRJNA1106903](https://www.ebi.ac.uk/ena/browser/view/PRJNA1106903).

```bash
mkdir -p ${PROJECT}/fastq
wget -nc -P ${PROJECT}/fastq ftp.sra.ebi.ac.uk/vol1/fastq/SRR288/058/SRR28867558/SRR28867558_1.fastq.gz
wget -nc -P ${PROJECT}/fastq ftp.sra.ebi.ac.uk/vol1/fastq/SRR288/058/SRR28867558/SRR28867558_2.fastq.gz
```

### Extract subset of 1M reads for profiling

```bash
mkdir -p ${PROJECT}/fastq/subset
FILES=(${PROJECT}/fastq/*.fastq.gz)
for FILE in ${FILES[@]}
do
    ID=$(basename ${FILE%.fast*})
    zcat -c ${FILE} | head --l 4000000 | gzip > ${PROJECT}/fastq/subset/${ID}.fastq.gz
done
```

### Generate barcode match profiles

This step requires barcode whitelists associated with the assay being used. Parse Bioscience customers can access the whitelists for the different assays by downloading their splitpipe pipeline. The whitelists are csv files in a barcodes directory (e.g. barcode_data_v1.csv). We only require the barcode sequence for scarecrow, so this needs cutting from the file (i.e. `cut -d',' -f2 barcode_data_v1.csv | sed '1d' > barcode_data_v1.txt`). Once the whitelists are generated, they can be defined as colon-delimited strings (`<barcode index>:<whitelist name>:<whitelist file>`) in a bash array for later use.

```bash
BARCODES=(BC1:n99_v5:${PROJECT}/barcode_whitelists/bc_data_n99_v5.txt
          BC2:v1:${PROJECT}/barcode_whitelists/bc_data_v1.txt
          BC3:v1:${PROJECT}/barcode_whitelists/bc_data_v1.txt)
```

We can now run `scarecrow seed` to process each barcode whitelist. The below example is for a SLURM HPC, but will work on a standard PC by omitting the `sbatch` line.

**check RAM usage before finalising documentation as it was previously 1G for 100K reads**

```bash
mkdir -p ${PROJECT}/barcode_profiles
THREADS=1
FASTQS=(${PROJECT}/fastq/subset/*.fastq.gz)
for BARCODE in ${BARCODES[@]}
do
    sbatch --ntasks ${THREADS} --mem 4G --time=01:00:00 -o seed.%j.out -e seed.%j.err \
        scarecrow seed --threads ${THREADS} \
            --fastqs ${FASTQS[@]} --strands pos neg \
            --barcodes ${BARCODE} \
            --out ${PROJECT}/barcode_profiles/barcodes.${BARCODE%%:*}.csv
done
```