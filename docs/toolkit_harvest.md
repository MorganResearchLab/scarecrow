<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[root](root.md)

## scarecrow harvest
The CSV file from each barcode whitelist are processed with `harvest` to identify barcode match peaks. The conserved sequences TSV file, which will be identical for any set of fastq files analysed with different barcode whitelists, can also be provided to mask peaks that fall within these regions. A minimum distance between barcodes can be provided, and the number of peaks (--barcode_count) to return per barcode whitelist. The results are recorded to a CSV file and a histogram of the peaks is plotted to a PNG file. *The CSV file should be edited if necessary to ensure that the barcode_whitelist values are unique*.

```bash
scarecrow harvest BC1_counts.csv BC2_counts.csv BC3_counts.csv \
    --barcode_count 1 \
    --min_distance 10 \
    --conserved BC1_conserved.tsv \
    --out barcode_positions.csv
```

The output file has the below format. `barcode_whitelist` indicates the whitelist from which the majority of barcodes at that barcode index can be found, `file_index` is the 0-based index of the FASTQ file (of FASTQ files passed to `harvest`) on which an exact match was identified; the associated FASTQ file for the index is listed under `file`. The `orientation` indicates the orientation of barcodes at that index on the sequencing read, `start` and `end` indicate the start and end positions of the barcode index, `read_count` indicates the number of reads with a barcode match at that location, while `read_fraction` presents the `read_count` as a fraction of the total number of reads analysed.

```bash
barcode_whitelist,file_index,file,orientation,start,end,read_count,read_fraction
BC2:v1,2,SRR28867558_3.fastq.gz,forward,11,18,8850,0.93
BC3:v1,2,SRR28867558_3.fastq.gz,forward,49,56,8327,0.87
BC1:n99_v5,2,SRR28867558_3.fastq.gz,forward,79,86,7447,0.8
```