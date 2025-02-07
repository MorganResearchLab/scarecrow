<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[Back to root](root.md)

# Example: 10X Flex

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../img/10Xflex_dark.svg">
  <img alt="Chromium 10X Flex library structure" src="../img/10Xflex_light.svg">
</picture>

Stay organised - create a folder for the project to keep things tidy.

```bash
PROJECT=./scarecrow/examples/10X3flex
mkdir -p ${PROJECT}
```

Download 10X Flex (FRP) data from [https://www.ebi.ac.uk/ena/browser/view/PRJNA1106903](https://www.ebi.ac.uk/ena/browser/view/PRJNA1106903).

```bash
mkdir -p ${PROJECT}/fastq
wget -nc -P ${PROJECT}/fastq ftp.sra.ebi.ac.uk/vol1/fastq/SRR288/063/SRR28867563/SRR28867563.fastq.gz
wget -nc -P ${PROJECT}/fastq ftp.sra.ebi.ac.uk/vol1/fastq/SRR288/062/SRR28867562/SRR28867562.fastq.gz
```

### 1. Generate barcode match profiles

