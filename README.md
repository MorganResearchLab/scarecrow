# scarecrow

<img style="float:right;width:200px;" src="./img/scarecrow.png" alt="scarecrow"/>

A library-agnostic toolkit for pre-processing combinatorial indexed single-cell RNA sequence data.

[Documentation](docs/root.md)

**scarecrow is undergoing substantial editing and may not behave as intended.**

### Todo

* Input validation
* Jitter does not currently apply to UMI or insert sequence
  - if UMI on same read and downstream then may need position updating before extraction
* Benchmark different assays (SPLiTseq, Parse, Scale) and methods (split-pipe, scarecrow, scaleRna)
  - barcode recovery
  - alignment
  - improved metrics (higher gene count, less variance, etc.)
  - DGE analysis with/without scarecrow


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
    time scarecrow seed --fastqs ${R1} ${R2} \
        -o ./WTv2/barcodes_${ID}_trie.csv --barcodes ${BARCODE} -n 0 -u 0 \
        --pickle ${WHITELIST}.${ID}.pkl.gz -k 2
done

# Harvest (set-based approach)
FILES=(./WTv2/barcodes.BC*.csv)
scarecrow harvest ${FILES[@]} --barcode_count 1 --min_distance 10 \
    --conserved ./WTv2/barcodes.BC1_conserved.tsv \
    --out ./WTv2/barcode_positions_set.csv

# Harvest (trie and kmer index approach)
FILES=(./WTv2/barcodes_BC*_trie.csv)
scarecrow harvest ${FILES[@]} --barcode_count 1 --min_distance 10 \
    --conserved ./WTv2/barcodes_BC1_trie_conserved.tsv \
    --out ./WTv2/barcode_positions_trie.csv

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
    --out ./cDNA_v2 --threads 1 --verbose &> debug.log

scarecrow weed --fastq P443A_index_10nt_1005_EKDL250000649-1A_22LJ3MLT4_L3_1.fq.gz \
    --sam cDNA_v2.sam \
    -i 1 \
    --out cDNA_v2_fix.sam \
    -m 1 \
    --barcodes BC3:P7:./BC3.txt &> debug.log
```



# Testing on laptop (10X3p)
```bash
R1=./10X3p/SRR28867562_3.1M.fastq.gz
R2=./10X3p/SRR28867562_4.1M.fastq.gz
BARCODE=(BC1:3M-Feb2018:./10X3p/3M-february-2018.txt)

# Generate custom trie and kmer index for use with scarecrow seed
time scarecrow encode --force_overwrite --barcodes ${BARCODE} --pickle -k 8
time scarecrow encode --barcodes ${BARCODE} --pickle -k 8

# Seed using trie
time scarecrow seed --fastqs ${R1} ${R2} \
    -o ./10X3p/barcodes_${BARCODE%:*:*}.csv \
    --barcodes ${BARCODE} -n 0 -u 0 \
    --pickle ./10X3p/3M-february-2018.txt.k8.pkl.gz

# Harvest (trie and kmer index approach)
FILES=(./10X3p/barcodes_BC1.csv)
scarecrow harvest ${FILES[@]} --barcode_count 1 --min_distance 10 \
    --conserved ./10X3p/barcodes_BC1_conserved.tsv \
    --out ./10X3p/barcode_positions.csv

# Reap (trie and kmer index approach)
BARCODE=(BC1:3M-Feb2018:./10X3p/3M-february-2018.txt.k8.pkl.gz)
time scarecrow reap --fastqs ${R1} ${R2} -j 0 -m 1 -q 10 \
    -p ./10X3p/barcode_positions.csv \
    --barcodes ${BARCODE} --extract 2:1-90 --umi 1:17-28 \
    --out ./10X3p/cDNA_k8 --threads 4 -b 10000

# Stats
scarecrow samstat --sam ./10X3p/cDNA_k8.sam
```
