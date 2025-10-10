<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[root](../README.md)

## scarecrow reap
The `reap` tool extracts a specified target sequence (e.g. cDNA on read 1 positions 1-64) and its associated quality values, together with cell barcodes and optionally the unique molecular index (UMI). The sequence data can be output to either an interleaved FASTQ file with `--out_fastq` or a SAM file with `--out_sam` (default), with the prefix `--out <prefix>`. If writing to SAM then the input FASTQ data should be trimmed beforehand to remove any adapter sequences. If output as a FASTQ file, a JSON file is also generated which contains the commandline parameters required to process the file with the kallisto-bustools workflow (kb count).

There are options to set `jitter`, which is the flanking distance to extend the barcode search from the start and end positions of the barcode peaks; `mismatch`, which is the maximum number of mismatches between an expected and observed barcode sequence. If `jitter` extends beyond the start of the sequence, for example if the barcode starts at position 1, then it will clip *n* end bases and insert *n* x `N` start bases when checking negative positions. The `N` bases count towards mismatches. Similarly, if `jitter` accounting for barcode length extends beyond the end of the sequence, then *n* start bases will be clipped and *n* x `N` bases will be appended to the end of the sequence. Any barcodes starting with a negative position will have Phred quality scores of 0 (!) assigned to bases in the matching barcode that are not present in the original sequence.

The default barcode matching strategy implemented in `scarerow` is a set-based method. When barcode whitelists contain thousands of barcodes the set-based method may not be feasible due to memory constraints. An alternative has also been implemented, which combines a trie-based method for exact matching and a k-mer indexing method for approximate-matching. This method requires the generation of a pickle file which contains the trie and k-mer index. The pickle file can generated when running the `seed` tool (see [Parse example](./example_evercode.md)) or it can be created with the `encode` [tool](./toolkit_encode.md).

The `base_quality` setting will mask bases as `N` if their Phred quality score is below the specified number, this is performed before barcode matching and can significantly reduce the number of matched barcodes if set too high. Where a barcode match is not identified the barcode is recorded as `NNNNNNN`. An additional flag, `--sift` can be used to skip writing any reads that have invalid barcodes (i.e. any corrected barcode containing an `N`).

Files containing barcode mismatch counts (`_mismatch_stats.csv`) and start positions (`_position_stats.csv`) counts are also generated.

This tool processes reads in batches and so generates temporary files in the process. These temporary files are automatically removed once the process is complete.

```bash
scarecrow reap --fastqs R1.fastq.gz R2.fastq.gz \
    --barcode_positions barcode_positions.csv \
    --barcodes BC1:v1:v1_whitelist.txt BC2:v1:v1_whitelist.txt BC3:n198:n198_whitelist.txt \
    --out_fastq \
    --out ./cDNA \
    --extract 1:1-64 \
    --umi 2:1-10 \
    --jitter 2 \
    --mismatches 2 \
    --base_quality 10 \
    --threads 16
```

### Example interleaved FASTQ format

The first read contains the barcodes and UMI, while the second read contains the `extract` sequence range. The barcodes and quality scores derive from the corrected barcodes. The sequence headers include the SAM read tags (see below) to facilitate the option to `recast` to SAM file format if required.

```bash
@SRR28867558.10002 1 VH01123:94:AACNK35M5:1:1101:31391:1625 CR:CTGGCATA_GAGCTGAA_CCTGTTGC:CY:CCCCCCC;_CCCCCCCC_CCCCCCCC:CB:CTGGCATA_GAGCTGAA_CCTGTTGC:XQ:CCCCCCC;_CCCCCCCC_CCCCCCCC:XP:11_49_79:XM:0_0_0:UR:CGCGGAGGTT:UY:CCCCCCCCCC/1
CTGGCATAGAGCTGAACCTGTTGCCGCGGAGGTT
+
CCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCC
@SRR28867558.10002 1 VH01123:94:AACNK35M5:1:1101:31391:1625 CR:CTGGCATA_GAGCTGAA_CCTGTTGC:CY:CCCCCCC;_CCCCCCCC_CCCCCCCC:CB:CTGGCATA_GAGCTGAA_CCTGTTGC:XQ:CCCCCCC;_CCCCCCCC_CCCCCCCC:XP:11_49_79:XM:0_0_0:UR:CGCGGAGGTT:UY:CCCCCCCCCC/2
GTTTCATATGTTGGCCAGGCTGGTCTCAAACTCCTGACCTCGTGAT
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCC
```

The accompanying JSON identifies the barcode and UMI positions for the interleaved fastq file, required for the custom -x string of `kb count`. In addition, the `kb count` string disables barcode whitelist checking (`-w NONE`) and indicates that the FASTQ file is interleaved (`--inleaved`).

```bash
{
    "description": "scarecrow",
    "barcodes": [
        {
            "range": "1:1-8",
            "whitelist": "v1_whitelist.txt"
        },
        {
            "range": "1:9-16",
            "whitelist": "v1_whitelist.txt"
        },
        {
            "range": "1:17-24",
            "whitelist": "n198_whitelist.txt"
        }
    ],
    "umi": [
        {
            "range": "1:25-34"
        }
    ],
    "kallisto-bustools": [
        {
            "kb count": "-i </path/to/transcriptome.idx> -g </path/to/transcripts_to_genes> -x 0,0,8,0,8,16,0,16,24:0,24,34:1,0,0 -w NONE --h5ad --inleaved -o <outdir> ./cDNA.fastq"
        }
    ]
}
```

### Example SAM format

The SAM output includes sequence tags for the barcodes and UMI, in addition to barcode start positions (XP) and mismatch counts (XM).

```bash
@HD     VN:1.6
@PG     PN:scarecrow    ID:reap VN:0.1.3        CL:--command recast --infile ./cDNA.fastq
SRR28867558.280001      4       *       0       0       *       *       0       0       CCGTCTCGCCCGCCGCGCCGGGGAGGTGGAGCACGAGCGCACGTGTTAGGACCCGAAAGATGGTGAACTAT CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC CR:Z:GACAGTGC_ACAGATTC_TATGTGTC CY:Z:CCCCCCCC_CCCCCCCC_CCCCCCCC CB:Z:GACAGTGC_ACAGATTC_TATGTGTC XQ:Z:CCCCCCCC_CCCCCCCC_CCCCCCCC XP:Z:11_49_79   XM:Z:0_0_0      UR:Z:TTCGTAGCAC UY:Z:CCCCCCCCCC
```

The sequence tags applied are listed below:

| TAG | Description |
| --- | ----------- |
| CR  | Uncorrected barcode sequence at expected position |
| CY  | Uncorrected barcode base qualities |
| CB  | Corrected barcode sequence accounting for jitter |
| XQ  | Corrected barcode base qualities (0[!] quality applied to bases in negative space) |
| XP  | Corrected barcode start position |
| XM  | Corrected barcode mismatch count |
| UR  | UMI sequence |
| UY  | UMI base qualities |

### Example mismatch_stats format

This reports the count of mismatches observed across reads. In the below example the data has 3 barcodes and has been processed with `--mismatches 2`. Reads with valid barcodes can therefore have from 0 to 6 (2 per barcode) mismatches. Reads with invalid barcodes can therefore contain between 1 and 3 invalid barcodes, indicated by the negative numbers.

```bash
mismatches,count
-3,3814125
-2,3590285
-1,14409173
0,125784335
1,11300683
2,6035706
3,1026087
4,823998
5,391577
6,298803
```

### Example position_stats format

This reports the count of of barcode start positions observed. In the below example the data has 3 barcodes and has been processed with `--jitter 1`. Invalid barcodes have a position of N, and this count should equal the absolute summed mismatch count of negative mismatch counts (i.e. `(3 x 3814125) + (2 x 3590285) + (1 x 14409173)`). We observe 3 peaks across the remaining positions, at 11, 49, and 79, corresponding with the expected start positions of the 3 barcodes. The asterisk (*) position indicates for those barcodes that more than one match for the same barcode sequence was found with the same number of mismatches and distance. In such instances, the barcode and number of mismatches are known however the exact position is not.

```bash
position,count
*,11950
10,3658309
11,154983669
12,1527298
48,9795569
49,146681602
50,2644153
78,13877626
79,133791532
80,2420490
N,33032118
```
