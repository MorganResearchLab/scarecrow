<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[root](root.md)

## scarecrow reap
The `reap` tool extracts a specified target sequence (e.g. cDNA on read 1 positions 1-64) and its associated quality values. This sequence data is written to a new fastq file. An optional sequencing read and range for the unique molecular identifier (UMI) can also be provided. The `jitter` value is the flanking distance to extend the barcode search from the start and end positions of the barcode peaks. `mismatch` is the maximum number of mismatches between an expected and observed barcode sequence. `base_quality` will mask bases as `N` if their Phred quality score is below the specified number, this is performed before barcode matching and can significantly reduce the number of matched barcodes if set too high. Where a barcode match is not identified the barcode is recorded as `null` in the sequence header.

```bash
# Reap target sequence from fastqs (TBC)
scarecrow reap --fastqs <paired_fastq_R1> <paired_fastq_R2> \
    --barcode_positions <barcode_positions.csv> \
    --barcodes BC1:v1:<v1_whitelist.txt> BC2:v1:<v1_whitelist.txt> BC3:n198:<n198_whitelist.txt> \
    --barcodes ${BARCODES[@]} \
    --out cDNA.fastq \
    --extract 1:1-64 \
    --umi 2:1-10 \
    --jitter 2 \
    --mismatches 2 \
    --base_quality 10 \
    --threads 4
```

Below is an example read written to the output fastq file. Here a `null` barcode has been recorded at barcode index `79` with a mismatch count of `-1`. This is due to the sequence at the barcode index not matching any of the barcodes on the whitelist after accounting for the number of permitted mismatches.

```bash
```

The sequence header includes a number of SAM tags that can be added to a SAM file using the `samtag` tool post-alignment:

| TAG | Description |
| --- | ----------- |
| CR  | Uncorrected barcode sequence |
| CY  | Barcode base qualities |
| CB  | Corrected barcode sequence |
| XP  | Barcode start position |
| XM  | Barcode pre-correction mismatch count |
| UR  | UMI sequence |
| UY  | UMI base qualities |
