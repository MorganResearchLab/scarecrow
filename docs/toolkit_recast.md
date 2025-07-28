<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[root](../README.md)

## scarecrow recast
The `recast` tool converts either a scareecrow SAM file to a scarecrow interleaved FASTQ with accompanying JSON file, or a scarecrow interleaved FASTQ file into a scarecrow SAM file. The first read in the scarecrow FASTQ contains the barcodes and UMI, whilst the second read contains the target (i.e. cDNA) sequence. The JSON file outlines the barcode positions and provides parameters for processing the FASTQ file with kallisto-bustools kb count. Barcode sequences, qualities, positions and mismatches, along with UMI sequence and barcodes, are recorded in the read header of the FASTQ file and the read tags of the SAM file. With regards to the duplication of some of this information in the FASTQ file, this is intentional to ensure cross-functionality of this tool and `scarecrow stats`. For more details on the output formats see [`scarecrow reap`](toolkit_reap.md).

```bash
scarecrow recast --in data.sam
scarecrow recast --in data.fastq
```
