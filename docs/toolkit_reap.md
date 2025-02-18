<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[root](root.md)

## scarecrow reap
The `reap` tool extracts a specified target sequence (e.g. cDNA on read 1 positions 1-64) and its associated quality values, together with cell barcodes and optionally the unique molecular index (UMI). The sequence data can be output to either a FASTQ file or a SAM file (default) with the prefix `--out`. If writing to SAM then it the input FASTQ data should be trimmed beforehand to remove any adapter sequences. `jitter` is the flanking distance to extend the barcode search from the start and end positions of the barcode peaks. `mismatch` is the maximum number of mismatches between an expected and observed barcode sequence. `base_quality` will mask bases as `N` if their Phred quality score is below the specified number, this is performed before barcode matching and can significantly reduce the number of matched barcodes if set too high. Where a barcode match is not identified the barcode is recorded as `NNNNNNN` in the sequence header. 

```bash
# Reap target sequence from fastqs (TBC)
scarecrow reap --fastqs <paired_fastq_R1> <paired_fastq_R2> \
    --barcode_positions <barcode_positions.csv> \
    --barcodes BC1:v1:<v1_whitelist.txt> BC2:v1:<v1_whitelist.txt> BC3:n198:<n198_whitelist.txt> \
    --barcodes ${BARCODES[@]} \
    --out ./cDNA \
    --extract 1:1-64 \
    --umi 2:1-10 \
    --jitter 2 \
    --mismatches 2 \
    --base_quality 10 \
    --threads 4
```

Below is an example read written in FASTQ format. Here a barcode has been corrected as `NNNNNNNN` at barcode position `51` with a mismatch count of `-1` indicating no barcode was found within the specified mismatch distance.

```bash
@LH00509:177:22W5HTLT3:1:1101:47416:1000 1:N:0:CAGATCAC+ATGTGAAG CR=GAGGCTGT_CATCAAGT_CCAGTTCA CY=IIIIIIII_IIIIIIII_IIIIIIII CB=NNNNNNNN_CATCAAGT_CCAGTTCA XP=51_31_11 XM=-1_0_0 UR=NGTTGTCTGT UY=#IIIIIIIII
AGCCGGCGGGAGCCNCGGGGAGAGTTCTCTTTTCTTTGTGAAGGGCAGGGCGCCCTGGAATGGG
+
IIIIIIIIIIIIII#IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

Below is the same read output in SAM format.

```bash
LH00509:177:22W5HTLT3:1:1101:47416:1000 4   *   0   255 *   *   0   0   AGCCGGCGGGAGCCNCGGGGAGAGTTCTCTTTTCTTTGTGAAGGGCAGGGCGCCCTGGAATGGG    IIIIIIIIIIIIII#IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    CR:Z:GAGGCTGT_CATCAAGT_CCAGTTCA CY:Z:IIIIIIII_IIIIIIII_IIIIIIII CB:Z:NNNNNNNN_CATCAAGT_CCAGTTCA XP:Z:51_31_11   XM:Z:-1_0_0 UR:Z:NGTTGTCTGT UY:Z:#IIIIIIIII
```

The sequence tags applied are listed below:

| TAG | Description |
| --- | ----------- |
| CR  | Uncorrected barcode sequence at expected position |
| CY  | Uncorrected barcode base qualities |
| CB  | Corrected barcode sequence accounting for jitter |
| XP  | Corrected barcode start position |
| XM  | Corrected barcode mismatch count |
| UR  | UMI sequence |
| UY  | UMI base qualities |
