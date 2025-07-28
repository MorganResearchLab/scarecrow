<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[root](../README.md)

## scarecrow stats
The `stats` tool reads either a SAM or FASTQ file and generates count data from the TAGS detailed below. It also counts the number of times each unique combination of uncorrected barcodes and corrected barcodes is observed.

```bash
scarecrow stats --in file.sam
```

| TAG | Description |
| --- | ----------- |
| CR  | Uncorrected barcode sequence at expected position |
| CB  | Corrected barcode sequence accounting for jitter |
| XP  | Corrected barcode start position |
| XM  | Corrected barcode mismatch count |
| UR  | UMI sequence |
