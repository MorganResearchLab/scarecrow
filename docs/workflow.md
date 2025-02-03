<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[Back to root](root.md)

# Workflow

The toolkit comprises several tools: `seed`, `harvest`, `reap`, `tally`, each of which feeds into the next. 

## scarecrow seed
Using a small subset of paired-end sequencing reads (e.g. 1M ), `seed` identifies barcode seeds within reads. It returns a CSV file of barcode alignments, a TSV file of nucleotide frequencies for each fastq file, and a TSV file of conserved sequence runs. This step should be repeated for each barcode whitelist. The barcode whitelist is a text file with one barcode sequence per line.

```bash
scarecrow seed --fastqs <fastq_R1> <fastq_R2> --strands pos neg --out <BC1_counts.csv> --barcodes BC1:v1:<v1_whitelist.txt> 
```

The `--out` file has the below format. The `read` column indicates the read on which a barcode match was found, `name` is the name retrieved from sequencing read, `barcode_whitelist` contains the first two elements of the string passed to `--barcdoes`, `barcode` is the sequence of the matched barcode, `orientation` is the orientation in which the sequence was found on the read, `start` and `end` are the positions of the alignment, `mismatches` is the number of mismatching bases between the fastq sequence and barcode.

```bash
read    name    barcode_whitelist       barcode orientation     start   end     mismatches
read1   SRR28867558.1   ('BC1', 'n99_v5')       CAATTTCC        forward 10      17      1
read1   SRR28867558.1   ('BC1', 'n99_v5')       AATCTTTC        reverse 19      26      1
read1   SRR28867558.1   ('BC1', 'n99_v5')       TACTGTCT        reverse 33      40      1
```

The `--out[-csv]_conserved.tsv` file has the below format. The `read` column indicates the read on which the conserved sequence was identified, `start` and `end` indicate the conserved sequence positions, `length` indicates its length, `sequence` is the conserved sequence, `median_frequency` is the median frequency of the sequence across the reads. Conserved sequences represent linkers, adapters, and other 'fixed' sequences within a read.

```bash
read    start   end     length  sequence        median_frequency
read2   31      48      18      TCGCATCGGCGTACGACT      0.8843
read2   57      78      22      ATCCACGTGCTTGAGACTGTGG  0.8721
```

The `--out[-csv]_read[1|2]_frequencies.tsv` file has the below format. The `read` column indicates the sequencing read analysed, the `position` indicates the position along the read, `A`, `C`, `G`, `T`, `N` indicate the frequencies of the observed nucleotides at that position across reads.

```bash
read    position        A       C       G       T       N
read1   1       0.4928  0.1722  0.2761  0.0590  0.0000
read1   2       0.4973  0.1531  0.1740  0.1747  0.0009
read1   3       0.2076  0.1803  0.4372  0.1749  0.0000
```

## scarecrow harvest
The CSV file from each barcode whitelist are processed with `harvest` to identify barcode match peaks. The conserved sequences TSV file, which will be identical for any set of fastq files analysed with different barcode whitelists, can also be provided to mask peaks that fall within these regions. A minimum distance between barcodes can be provided, and the number of peaks (--barcode_count) to return per barcode whitelist. The results are recorded to a CSV file and a histogram of the peaks is plotted to a PNG file. *The CSV file should be edited if necessary to ensure that the barcode_whitelist values are unique*.

```bash
scarecrow harvest <BC1_counts.csv> <BC2_counts.csv> <BC3_counts.csv> \
    --barcode_count 3 --min_distance 10 --conserved <BC1_conserved.tsv> \
    --out <barcode_positions.csv>
```

The output file has the below format. `barcode_whitelist` indicates the whitelist from which the majority of barcodes at that barcode index can be found, `read` indicates the sequencing read on which the barcode index is located, `orientation` indicates the orientation of barcodes at that index on the sequencing read, `start` and `end` indicate the start and end positions of the barcode index, `read_count` indicates the number of reads with a barcode match at that location, `read_fraction` presents the `read_count` as a fraction of the total number of reads analysed, `barcode_diversity` is not currently used but indicates the diversity of barcode sequences at that index.

```bash
barcode_whitelist,read,orientation,start,end,read_count,read_fraction,barcode_diversity
"('BC2', 'v1')",read2,forward,11,18,911491,0.93,0.0001
"('BC3', 'v1')",read2,forward,49,56,861490,0.88,0.0001
"('BC1', 'n99_v5')",read2,forward,79,86,788817,0.8,0.0001
```

## scarecrow reap
The `reap` tool extracts a specified target sequence (e.g. cDNA) and its associated quality values. This sequence data is written to a new fastq file, and the combination of cell barcodes identified near the positions previously identifed by `harvest` are appened to the sequence header. An optimal sequence range for the unique molecular identifier (UMI) can also be provided, which will also append the UMI sequence to the header. The `jitter` value is the flanking distance to extend the barcode search from the start and end positions of the barcode peaks. `mismatch` is the maximum number of mismatches between an expected and observed barcode sequence. `base_quality` will mask bases as `N` if their Phred quality score is below the specified number, this is performed before barcode matching and can significantly reduce the number of matched barcodes if set too high. Where a barcode match is not identified the barcode is recorded as `null` in the sequence header.

```bash
# Reap target sequence from fastqs (TBC)
scarecrow reap --fastqs <paired_fastq_R1> <paired_fastq_R2> --barcode_positions <barcode_positions.csv> \
    --barcodes BC1:v1:<v1_whitelist.txt> BC2:v1:<v1_whitelist.txt> BC3:n198:<n198_whitelist.txt> \
    -j 2 -m 2 -q 30 --barcodes ${BARCODES[@]} --extract 1:1-64 --umi 2:1-10 --out ./cDNA.fq --threads 4
    --jitter 2 --mismatches 2 --base_quality 10 --extract 1:1-64 --out cDNA.fastq
```

Below is an example read written to the output fastq file. Here a `null` barcode has been recorded at barcode index `79` with a mismatch count of `-1`. This is due to the sequence at the barcode index not matching any of the barcodes on the whitelist after accounting for the number of permitted mismatches.

```bash
```

| TAG | Description |
| --- | ----------- |
| CR  | Uncorrected barcode sequence |
| CY  | Barcode base qualities |
| CB  | Corrected barcode sequence |
| XP  | Barcode start position |
| XM  | Barcode pre-correction mismatch count |
| UR  | UMI sequence |
| UY  | UMI base qualities |

## scarecrow tally
The `tally` tool processes the fastq file generated by `reap` to record some metrics on the barcodes identified. The outputs include: `.barcodes.csv`, `.barcode.combinations.csv`, `.mismatches.csv`, `.positions`.

The `.barcodes.csv` file records at each barcode index the barcode sequence identified and the number of times it is observed.

```bash
Index,Barcode,Count
1,TCATTCGC,52
1,GTCTGCTC,818
1,GTCTCTGC,826
```

The `.barcode.combinations.csv` file records for each combination of barcodes the number of times the combination is observed.

```bash
BarcodeCombination,Count
CTGACTTC_CCTAATCC_GATAGACA,1
CATGTCTC_GTGTTCTA_GCCACATA,1
GCTATCAT_CCGAAGTA_CCGAAGTA,1
```

The `.mismatches.csv` file records the count of barcodes (not reads) with a given number of barcode mismatches.

```bash
BarcodeMismatches,Count
0,2432376
1,87246
```

The `.positions.csv` file records for each barcode index the numnber of barcodes found to start at a given position. For instance, if the barcode is expected at position 78, and jitter = 2, then the start position could be in the range 76-80.

```bash
Index,Position,Count
1,77,18226
1,78,74607
1,79,907167
```

[Back to root](root.md)