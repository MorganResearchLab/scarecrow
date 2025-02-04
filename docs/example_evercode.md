<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[Back to root](root.md)

# Example: Parse Evercode WTv2

#### Download data

Download data from [https://www.ebi.ac.uk/ena/browser/view/PRJNA1106903](https://www.ebi.ac.uk/ena/browser/view/PRJNA1106903).

```bash
mkdir -p ./fastq
wget -nc -P ./fastq ftp.sra.ebi.ac.uk/vol1/fastq/SRR288/058/SRR28867558/SRR28867558_1.fastq.gz
wget -nc -P ./fastq ftp.sra.ebi.ac.uk/vol1/fastq/SRR288/058/SRR28867558/SRR28867558_2.fastq.gz
```