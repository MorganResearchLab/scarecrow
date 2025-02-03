<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[root](root.md) | [toolkit](root.md)

## scarecrow seed
Using a small subset of paired-end sequencing reads (e.g. 1M ), `seed` identifies barcode seeds within reads. It returns a CSV file of barcode alignments, a TSV file of nucleotide frequencies for each fastq file, and a TSV file of conserved sequence runs. This step should be repeated for each barcode whitelist. The barcode whitelist is a text file with one barcode sequence per line.

```bash
scarecrow seed --fastqs <paired_fastq_R1> <paired_fastq_R2> \
    --strands pos neg \
    --out <BC1_counts.csv> \
    --barcodes BC1:v1:<v1_whitelist.txt> 
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