# scarecrow

![scarecrow](img/scarecrow.png)

A toolkit for preprocessing single cell sequencing data.

[Documentation](docs/root.md)

### Todo

* Error handling implementing to capture missing or incorrect parameters, and unexpected file content
* Peaks in between barcodes need further investigation
* Plot generated by harvest currently will not handle > 1 barcode peak per whitelist (doesn't affect CSV output)
* Benchmark different assays (SPLiTseq, Parse, 10X) and methods (split-pipe, scarecrow, UMI tools)
*   - barcode recovery
*   - alignment (STAR and kallisto)
* Test alignment with kallisto and STAR
*    - may need to alter sequence header formatting depending on what is retained in BAM file



### Parse Evercode WTv3 Example
```bash
OUTDIR=/uoa/home/s14dw4/scarecrow_test/MorganLab

# seed
R1=/uoa/scratch/shared/Morgan_Lab/Parse/LV6000778329-SO8229-University-of-Aberdeen-Morgan-Lab-Lib1-2024-12-10_S1_L001_R1_001.fastq.gz
R2=/uoa/scratch/shared/Morgan_Lab/Parse/LV6000778329-SO8229-University-of-Aberdeen-Morgan-Lab-Lib1-2024-12-10_S1_L001_R2_001.fastq.gz

# barcodes
BARCODES=(
  "BC1:R1_v3:${OUTDIR}/barcodes/bc_data_n123_R1_v3_5.csv"
  "BC2:v1:${OUTDIR}/barcodes/bc_data_v1.csv"
  "BC3:R3_v3:${OUTDIR}/barcodes/bc_data_R3_v3.csv"
  )

# Exract 100K reads
zcat -c ${R1} | head --l 400000 > ${OUTDIR}/fastq/100K_1.fastq 
zcat -c ${R2} | head --l 400000 > ${OUTDIR}/fastq/100K_2.fastq 

# Identify barcode seeds
for BARCODE in ${BARCODES[@]}
do
    scarecrow seed --fastqs ${R1} ${R2} --strands pos neg \
      -o ${OUTDIR}/results/barcodes_${BARCODE%:*:*}.csv --barcodes ${BARCODE}
done    

# Harvest barcode peaks
FILES=(${OUTDIR}/results/barcodes_*csv)
scarecrow harvest ${FILES[@]} --barcode_count 3 --min_distance 11 \
    --out ${OUTDIR}/barcode_positions.csv

# Reap cDNA
scarecrow reap --fastqs ${R1} ${R2} -p ${OUTDIR}/barcode_positions.csv \
    -j 5 -m 1 --barcodes ${BARCODES[@]} --read1 0-64 --out ./cDNA.fq.gz  
```






# Testing on laptop
```bash
R1=100K_1.fastq
R2=100K_2.fastq
BARCODES=(BC1:R1_v3:/Users/s14dw4/Documents/scarecrow_test/barcodes/bc_data_n123_R1_v3_5.barcodes
          BC2:v1:/Users/s14dw4/Documents/scarecrow_test/barcodes/bc_data_v1.barcodes
          BC3:R3_v3:/Users/s14dw4/Documents/scarecrow_test/barcodes/bc_data_R3_v3.barcodes)
for BARCODE in ${BARCODES[@]}
do
    scarecrow seed --fastqs ${R1} ${R2} --strands pos neg \
      -o ./results/barcodes_${BARCODE%:*:*}.csv --barcodes ${BARCODE}
done

FILES=(./results/barcodes_BC*csv)
scarecrow harvest ${FILES[@]} --barcode_count 3 --min_distance 11 \
    --conserved ./results/barcodes_BC1_conserved.tsv --out barcode_positions.csv

time scarecrow reap --fastqs ${R1} ${R2} -p ./barcode_positions.csv --barcode_reverse_order \
    -j 2 -m 2 -q 30 --barcodes ${BARCODES[@]} --extract 1:1-64 --umi 2:1-10 --out ./cDNA.fq --threads 4

scarecrow tally -f ./cDNA.fq -m 2

```




# Test issue with barcode not being matched (this is due to low quality scores)
# @SRR28867558.1 2 VH01123:94:AACNK35M5:1:1101:65343:1000/3
# TGTTATAACC[ATCATTCC]GTGGCCGATGTTTCGCATCGGCGTACGACT[AAGGTACA]ATCCACGTGCTTGAGACTGTGG[TTACCTGC]
# @SRR28867558.1 1 VH01123:94:AACNK35M5:1:1101:65343:1000/2 barcodes=TTACCTGC_AAGGTACA_null positions=79_49_11 mismatches=0_0_-1 UMI=TGTTATAACC
```bash
OUT=/Users/s14dw4/Documents/scarecrow_test/WTv2
FASTQS=(${OUT}/*fastq.gz)
for FASTQ in ${FASTQS[@]}
do
    gunzip -c ${FASTQ} | head --l 400 > ${FASTQ%.gz}
done
POS=${OUT}/barcode_positions.csv
BARCODES=(BC1:n99_v5:${OUT}/bc_data_n99_v5.txt
          BC2:v1:${OUT}/bc_data_v1.txt
          BC3:v1:${OUT}/bc_data_v1.txt)
time scarecrow reap --fastqs ${FASTQS[@]%.gz} -p ${POS} --barcode_reverse_order \
    -j 2 -m 2 -q 20 --barcodes ${BARCODES[@]} --extract 1:1-64 --umi 2:1-10 \
    --out ${OUT}/cDNA.fq --threads 4 
```





Details on RG tags:
https://github.com/samtools/hts-specs
https://github.com/samtools/hts-specs/blob/master/SAMtags.pdf
https://www.biostars.org/p/9593008/

```python
import pysam

def add_tags_to_sam(input_sam, output_sam):
    with pysam.AlignmentFile(input_sam, "r") as infile, pysam.AlignmentFile(output_sam, "w", header=infile.header) as outfile:
        for read in infile:
            # Extract the sequence header
            fastq_header = read.query_name
            
            # Parse the barcodes and UMI from the FASTQ header
            # Example header: "@LH00509:177:22W5HTLT3:1:1101:46251:1000 1:N:0:CAGATCAC+ATGTGAAG barcodes=TCTGATCC_GAACAGGC_ATCCTGTA positions=51_31_11 mismatches=0_0_0 UMI=NGAACTGAGT"
            try:
                details = fastq_header.split(" ")
                attributes = {item.split("=")[0]: item.split("=")[1] for item in details if "=" in item}
                
                # Add the CB and ZU tags if available
                if "barcodes" in attributes:
                    read.set_tag("CB", attributes["barcodes"], value_type="Z")
                if "UMI" in attributes:
                    read.set_tag("ZU", attributes["UMI"], value_type="Z")
            except Exception as e:
                print(f"Failed to parse header for read {fastq_header}: {e}")

            # Write the updated read to the output SAM file
            outfile.write(read)

# Input and output SAM file paths
input_sam = "input.sam"
output_sam = "output.sam"

add_tags_to_sam(input_sam, output_sam)
```

More efficient approach, which may require a lot of RAM

```python
import pysam

def preprocess_fastq_headers(fastq_file):
    """Preprocess FASTQ headers to create a mapping of read names to tags."""
    read_tags = {}
    with open(fastq_file, "r") as fq:
        for line in fq:
            if line.startswith("@"):
                # Extract the read name and tags
                fastq_header = line.strip()
                read_name = fastq_header.split(" ")[0][1:]  # Remove "@" and split to get the read name
                attributes = {item.split("=")[0]: item.split("=")[1] for item in fastq_header.split() if "=" in item}
                
                # Store barcodes and UMI if they exist
                barcodes = attributes.get("barcodes", None)
                umi = attributes.get("UMI", None)
                read_tags[read_name] = {"CB": barcodes, "ZU": umi}
    return read_tags

def add_tags_to_sam(input_sam, output_sam, read_tags):
    """Add tags to SAM file based on the preprocessed FASTQ headers."""
    with pysam.AlignmentFile(input_sam, "r") as infile, pysam.AlignmentFile(output_sam, "w", header=infile.header) as outfile:
        for read in infile:
            if read.query_name in read_tags:
                tags = read_tags[read.query_name]
                if tags["CB"]:
                    read.set_tag("CB", tags["CB"], value_type="Z")
                if tags["ZU"]:
                    read.set_tag("ZU", tags["ZU"], value_type="Z")
            outfile.write(read)

# Paths to the input files
fastq_file = "input.fastq"
input_sam = "input.sam"
output_sam = "output.sam"

# Preprocess FASTQ headers and add tags
read_tags = preprocess_fastq_headers(fastq_file)
add_tags_to_sam(input_sam, output_sam, read_tags)
```
