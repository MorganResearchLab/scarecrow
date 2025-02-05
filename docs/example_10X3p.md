<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[Back to root](root.md)

# Example: 10X 3'

<img style="float:right;width:100%;" src="../img/10X3p.svg" alt="scarecrow"/>

Stay organised - create a folder for the project to keep things tidy.

```bash
PROJECT=./scarecrow/examples/10X3p
mkdir -p ${PROJECT}
```

Download 10X 3' data from [https://www.ebi.ac.uk/ena/browser/view/PRJNA1106903](https://www.ebi.ac.uk/ena/browser/view/PRJNA1106903).

```bash
mkdir -p ${PROJECT}/fastq
wget -nc -P ${PROJECT}/fastq ftp.sra.ebi.ac.uk/vol1/fastq/SRR288/063/SRR28867563/SRR28867563.fastq.gz
wget -nc -P ${PROJECT}/fastq ftp.sra.ebi.ac.uk/vol1/fastq/SRR288/062/SRR28867562/SRR28867562.fastq.gz
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

Chromium barcode whitelists for different chemistry versions are available at [https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html).