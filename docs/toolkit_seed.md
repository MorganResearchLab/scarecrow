<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[root](../README.md)

## scarecrow seed
This tool identifies barcode seed positions and potential linker sequences within reads. It returns a CSV file of barcode alignments, a TSV file of nucleotide frequencies for each fastq file, and a TSV file of conserved sequence runs. This step should be repeated for each barcode whitelist. The barcode whitelist is a text file with one barcode sequence per line and no header. The `--num_reads` option specifies the number of reads to run the analysis on, and defaults to 100000. Using all reads (`--num reads 0`) is neither necessary nor recommended, the tool does not run multiprocessing so will take considerable time to process millions of reads. Subsets of reads are sampled randomly, the random seed can be set by `--random_seed <int>` and defaults to 1234.

```bash
scarecrow seed --fastqs R1.fastq.gz R2.fastq.gz \
    --num_reads 100000 \
    --out BC1_counts.csv \
    --barcodes BC1:v1:v1_whitelist.txt
```

The `--out` file has the below format. The `file_index` is the 0-based index of the FASTQ file (of files passed via `--fastq`) on which an exact match was identified; the associated FASTQ file for the index is listed under `file`. The `read_name` column indicates the read name of the sequencing read, and `seqlen` indicates the length of the read. The `barcode_whitelist` contains the first two elements of the string passed to `--barcdoes`, `barcode` is the sequence of the matched barcode, `orientation` is the orientation in which the sequence was found on the read, `start` and `end` are the positions of the alignment, `mismatches` is the number of mismatching bases between the fastq sequence and barcode.

```bash
file_index      file    read_name       seqlen  barcode_whitelist       barcode orientation     start   end     mismatches
1       SRR28867558_2.fastq.gz  SRR28867558.5   74      BC1:n99_v5      CACTTTCA        reverse 8       15      0
1       SRR28867558_2.fastq.gz  SRR28867558.5   74      BC1:n99_v5      GTGCTTGA        reverse 24      31      0
1       SRR28867558_2.fastq.gz  SRR28867558.5   74      BC1:n99_v5      GTGCTTGA        reverse 50      57      0
2       SRR28867558_3.fastq.gz  SRR28867558.18  86      BC1:n99_v5      GTGCTTGA        forward 63      70      0
2       SRR28867558_3.fastq.gz  SRR28867558.18  86      BC1:n99_v5      TGTGTATG        forward 79      86      0
```

The `--out[-csv]_conserved.tsv` file has the below format. The `file_index` column indicates the 0-based FASTQ file index from which the conserved sequence was identified, `start` and `end` indicate the conserved sequence positions, `length` indicates its length, `sequence` is the conserved sequence, `median_frequency` is the median frequency of the sequence across the reads. Conserved sequences represent linkers, adapters, and other 'fixed' sequences within a read.

```bash
file_index      start   end     length  sequence        median_frequency
2       31      48      18      TCGCATCGGCGTACGACT      0.8837
2       57      78      22      ATCCACGTGCTTGAGACTGTGG  0.8698
```

The `--out[-csv]_read[1|2]_frequencies.tsv` file has the below format. The `read` column indicates the sequencing read analysed, the `position` indicates the position along the read, `A`, `C`, `G`, `T`, `N` indicate the frequencies of the observed nucleotides at that position across reads.

```bash
read    position        A       C       G       T       N
read1   1       0.4928  0.1722  0.2761  0.0590  0.0000
read1   2       0.4973  0.1531  0.1740  0.1747  0.0009
read1   3       0.2076  0.1803  0.4372  0.1749  0.0000
```
