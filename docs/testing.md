# Testing on laptop (WTv2)
```bash
R1=./WTv2/100K_R1.fastq
R2=./WTv2/100K_R2.fastq
BARCODES=(BC1:n99_v5:./WTv2/bc_data_n99_v5.txt
          BC2:v1:./WTv2/bc_data_v1.txt
          BC3:v1:./WTv2/bc_data_v1.txt)

# Seed
for BARCODE in ${BARCODES[@]}
do
    ID=${BARCODE%:*:*}
    WHITELIST=${BARCODE#*:*:}
    echo ${ID}
    time scarecrow seed --fastqs ${R1} ${R2} \
        -o ./WTv2/barcodes_${ID}_set.csv --barcodes ${BARCODE} -n 0 -u 0
done

# Harvest (set-based approach)
FILES=(./WTv2/barcodes.BC*.csv)
scarecrow harvest ${FILES[@]} --barcode_count 1 --min_distance 10 \
    --conserved ./WTv2/barcodes.BC1_conserved.tsv \
    --out ./WTv2/barcode_positions_set.csv

# Reap (set-based approach)
BARCODES=(BC1:n99_v5:./WTv2/bc_data_n99_v5.txt
          BC2:v1:./WTv2/bc_data_v1.txt
          BC3:v1:./WTv2/bc_data_v1.txt)
time scarecrow reap --fastqs ${R1} ${R2} -j 1 -m 2 -q 10 \
    -p ./WTv2/barcode_positions_set.csv \
    --barcodes ${BARCODES[@]} --extract 1:1-74 --umi 2:1-10 \
    --out ./WTv2/cDNA --threads 2 --out_sam

# Reap (trie and kmer index approach)
BARCODES=(BC1:n99_v5:./WTv2/bc_data_n99_v5.txt.BC1.pkl.gz
          BC2:v1:./WTv2/bc_data_v1.txt.BC2.pkl.gz
          BC3:v1:./WTv2/bc_data_v1.txt.BC3.pkl.gz)
time scarecrow reap --fastqs ${R1} ${R2} -j 1 -m 2 -q 10 \
    -p ./WTv2/barcode_positions_trie.csv \
    --barcodes ${BARCODES[@]} --extract 1:1-74 --umi 2:1-10 \
    --out ./WTv2/cDNA_trie --threads 2

scarecrow samstat --sam ./WTv2/cDNA_set.sam
scarecrow samstat --sam ./WTv2/cDNA_trie.sam


mamba activate kallisto
pip install kallisto
pip install kb-python
kb ref -i ref/transcriptome.idx -g ref/transcripts_to_genes.txt -d human

mkdir ./WTv2/kb

# Whitelist not working, kallisto not getting correct barcode length 
#cat ./WTv2/bc_data_n99_v5.txt ./WTv2/bc_data_v1.txt ./WTv2/bc_data_v1.txt > ./WTv2/WTv2_cmb.txt
kb count -i ./ref/transcriptome.idx -g ./ref/transcripts_to_genes.txt \
    -x 0,0,7,0,8,15,0,16,23:0,24,33:1,0,0 -w NONE --h5ad \
    --overwrite --inleaved -o ./WTv2/kb ./WTv2/cDNA_set.fastq


# needs to run on an aligned file
mkdir ./WTv2/umi_tools
umi_tools dedup --stdin ./WTv2/cDNA_set.sam --output-stats=./WTv2/umi_tools/cDNA \
    --extract-umi-method=tag --umi-tag=UR --cell-tag=CB


```




# Testing on laptop (Scale)
```bash
FASTQS=(./Scale/*fastq.gz)
BARCODES=(BC1:lig:./Scale/3lvlRNA_lig.txt
          BC2:rt:./Scale/3lvlRNA_rt.txt
          BC3:pcr:./Scale/3lvlRNA_pcr.txt)

# Seed
for BARCODE in ${BARCODES[@]}
do
    ID=${BARCODE%:*:*}
    WHITELIST=${BARCODE#*:*:}
    echo ${ID}
    time scarecrow seed --fastqs ${FASTQS[@]} \
        -o ./Scale/barcodes_${ID}.csv --barcodes ${BARCODE} -n 0 -u 0
done

# Harvest (set-based approach)
FILES=(./Scale/barcodes_BC*.csv)
scarecrow harvest ${FILES[@]} --barcode_count 1 --min_distance 10 \
    --conserved ./Scale/barcodes_BC1_conserved.tsv \
    --out ./Scale/barcode_positions.csv


# scarecrow reap
BQ=10
JITTER=1
MISMATCH=1
scarecrow reap --threads 2 --batch_size 50000 \
                --fastqs ${FASTQS[@]} \
                -p ./Scale/barcode_positions.csv \
                --barcode_reverse_order \
                --barcodes ${BARCODES[@]} \
                --extract 4:1-76 --umi 3:16-23 \
                --jitter ${JITTER} --mismatch ${MISMATCH} --base_quality ${BQ} \
                --out ./Scale/cDNA




# Test fastq output (qual v seq len)
BQ=10
JITTER=2
MISMATCH=2
FASTQS=(./test/*.fastq)
BARCODES=(BC1:n141_R1_v3_6:./barcodes/Parse/bc_data_n141_R1_v3_6.txt
          BC2:v1:./barcodes/Parse/bc_data_v1.txt
          BC3:R3_v3:./barcodes/Parse/bc_data_R3_v3.txt)
scarecrow reap --threads 1 --batch_size 50000 \
    --fastqs ${FASTQS[@]} \
    -p ./barcode_profiles/Parse-Morgan/barcode_positions.csv \
    --barcodes ${BARCODES[@]} \
    --extract 1:1-64 --umi 2:1-10 \
    --jitter ${JITTER} --mismatch ${MISMATCH} --base_quality ${BQ} \
    --out ./test/cDNA \
    --out_fastq

```



# Testing on laptop (split-seq)
```bash
cd ~/Documents/split-seq
R1=./r1.fastq
R2=./r2.fastq
BARCODES=(BC1:3lvl:./BC1.txt
          BC2:3lvl_lig:./BC2.txt
          BC3:P7:./BC3.txt)

# Seed
for BARCODE in ${BARCODES[@]}
do
    ID=${BARCODE%:*:*}
    WHITELIST=${BARCODE#*:*:}
    echo ${ID}
    time scarecrow seed --fastqs ${R1} ${R2} \
        -o ./barcodes_${ID}.csv --barcodes ${BARCODE} -n 0 -u 0
done

FILES=(./barcodes_BC*.csv)
scarecrow harvest ${FILES[@]} --barcode_count 1 --min_distance 10 \
    --conserved ./barcodes_BC1_conserved.tsv \
    --out ./barcode_positions.csv

scarecrow reap --fastqs ${R1} ${R2} -j 2 -m 3 -q 10 \
    -p ./barcode_positions.csv \
    --barcodes ${BARCODES[@]} --extract 2:11-150 --umi 2:1-10 --base_quality 10 \
    --out ./cDNA_v2 --threads 1 

scarecrow weed --fastq P443A_index_10nt_1005_EKDL250000649-1A_22LJ3MLT4_L3_1.fq.gz \
    --sam cDNA_v2.sam \
    -i 1 \
    --out cDNA_v2_fix.sam \
    -m 1 \
    --barcodes BC3:P7:./BC3.txt &> debug.log
```
