<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[root](../README.md)

## scarecrow rake
The `rake` tool injects a scarecrow corrected barcode from the CB tag of a scarecrow SAM file into a source FASTQ file at a specified position. This is intended for use on 10X sequencing data. The principle being that Cell Ranger is first run on the FASTQ files. Reads that were not assigned a barcode are then extracted using the `seqkit grep` [tool](https://bioinf.shenwei.me/seqkit/) and are processed with scarecrow (seed -> harvest -> reap) to match barcodes accounting for both jitter and mismatches. The SAM file returned by scarecrow is then provided to `rake`, along with the 10X FASTQ containing the barcode and UMI sequence (e.g. R1.fastq.gz). The tool will then retrieve reads with valid barcodes from the SAM file, and replace the sequence starting at position x in the FASTQ file with the CB sequence. The tool defaults to position 1 with barcode length 16, as per the 10X library. Cell Ranger is then re-run using the new FASTQ file in place of the original - note you will need to ensure the file is renamed to suit the original format.

```bash
scarecrow rake --sam scarecrow.sam --fastq 10X_R1.fastq.gz --out 10X_R1.rake.fastq.gz --position 1 --length 16
```

The read name, CR (uncorrected barcode) and CB (corrected barcode) tags can be extracted from the CellRanger possorted_genome BAM file using samtools and awk as follows:

```bash
samtools view --threads ${THREADS} ${SAM} | awk '{
    read_name = $1;
    tags = "";
    for (i = 12; i <= NF; i++) {
        if ($i ~ /^CB:Z:/ || $i ~ /^CR:Z:/) {
            tags = tags "\t" $i;
        }
    }
    print read_name tags;
}' > ${SAM}.tags
```

Which will return a `.tags` file similar to the below. Reads that have no CB tag are those that were not assigned a barcode by CellRanger.

```bash
SRR28867563.25150069	CR:Z:ATCGTCCTCTAACGCA	CB:Z:ATCGTCCTCTAACGCA-1
SRR28867563.57692158	CR:Z:TGAATGCCACGTATAC	CB:Z:TGAATGCCACGTATAC-1
SRR28867563.110363081	CR:Z:GATCATGTCCAGCTCT	CB:Z:GATCATGTCCAGCTCT-1
SRR28867563.108926426	CR:Z:GTACAGTGTGGGTCAA	CB:Z:GTACAGTGTGGGTCAA-1
SRR28867563.112330373	CR:Z:CATACCCCACAGCCTG	CB:Z:CATACCCCACAGCCTG-1
SRR28867563.188967795	CR:Z:CCTGCATTCGCAATGT	CB:Z:CCTGCATTCGCAATGT-1
SRR28867563.94538128	CR:Z:TTGCTGCCTCTTAGCC
```

The reads containing missing and CellRanger-corrected barcodes can be written to separate files as follows:

```bash
awk -F'\t' -v missing_cb_file=missing_CB.txt \
    -v cb_mismatch_file=CB_not_CR.txt \
    '
    {
        # Extract CR sequence
        cr = $2
        sub(/^CR:Z:/, "", cr)
        if (NF < 3) {
            # Missing CB
            missing_cb++
            print $1 >> missing_cb_file
        } else {
            # Extract CB sequence
            cb = $3
            sub(/^CB:Z:/, "", cb)
            sub(/-[0-9]+$/, "", cb)

            if (cb != cr) {
                cb_mismatch++
                print $1 >> cb_mismatch_file
            }
        }
    }
    END {
        print "Missing CB:", missing_cb
        print "CB != CR:", cb_mismatch
    }
    ' "${SAM}.tags"
```

The reads missing barcodes can then be extracted from the original FASTQ files using seqkit as follows:

```bash
FASTQS=(./fastq/*fastq.gz)
mkdir -p ./fastq/missing
for FASTQ in ${FASTQS[@]}
do
    ID=$(basename ${FASTQ})
    seqkit grep \
        -f missing_CB.txt \
        -o ./fastq/missing/${ID} \
        ${FASTQ}
done
```

These FASTQ files containing the subset of reads can then be processed with scarecrow (seed, harvest, reap) to output a SAM file containing corrected barcodes accounting for jitter and mismatches. The resulting SAM file is used with rake, alongside the original FASTQ file containing the barcode, to generate a corrected FASTQ file for re-use with CellRanger.
