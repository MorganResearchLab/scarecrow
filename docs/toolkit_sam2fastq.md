<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[root](root.md)

## scarecrow sam2fastq
The `sam2fastq` tool converts a SAM file to a scarecrow interleaved FASTQ with accompanying JSON file. The SAM file need not be aligned, and need not contain @SQ header lines. The first read in the FASTQ contains the barcodes and UMI, whilst the second read contains the read sequence. The JSON file outlines the barcode positions and provides parameters for processing the FASTQ file with kallisto-bustools kb count. For more details see [`scarecrow reap`](toolkit_reap.md).

```bash
scarecrow sam2fastq --sam aligned.bam
```
