<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[Back to root](root.md)

# Example: 10X 5'

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../img/10X5p_dark.svg">
  <img alt="Chromium 10X 5' library structure" src="../img/10X5p_light.svg">
</picture>

Stay organised - create a folder for the project to keep things tidy.

```bash
PROJECT=./scarecrow/examples/10X5p
mkdir -p ${PROJECT}
```

Download 10X 5' data from [https://www.ebi.ac.uk/ena/browser/view/PRJNA1106903](https://www.ebi.ac.uk/ena/browser/view/PRJNA1106903).

```bash
mkdir -p ${PROJECT}/fastq
wget -nc -P ${PROJECT}/fastq ftp.sra.ebi.ac.uk/vol1/fastq/SRR288/063/SRR28867563/SRR28867563.fastq.gz
wget -nc -P ${PROJECT}/fastq ftp.sra.ebi.ac.uk/vol1/fastq/SRR288/062/SRR28867562/SRR28867562.fastq.gz
```

### 1. Generate barcode match profiles

Chromium barcode whitelists for different chemistry versions are available at [https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium5.html](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium5.html).