<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[root](../README.md)

## scarecrow pickle
This tool generates an Aho-Corasick trie and either a seed index (default) or k-mer index from a whitelist of barcodes. The seed index applies the pigeonhole principle to determine the number of seeds based on the number of mismatches allowed, whilst the k-mer index chunks the barcode into k = L // (m+1) overlapping k-mers. The resulting trie and index is written to a file with the input whtelist filename as the prefix with the additional suffix of `.${index_type}.m{mismatches}.pkl.gz`.

```bash
scarecrow pickle --barcodes BC1:v1:v1_whitelist.txt --index kmer --mismatches 2
```
