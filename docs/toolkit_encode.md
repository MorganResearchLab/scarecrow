<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[root](../README.md)

## scarecrow encode
This tool generates an Aho-Corasick trie and either a seed index (default) or k-mer index from a whitelist of barcodes. Whichever index is used, the pidgeonhole principle is applied to determine the number of chunks (seeds or k-mer size) required for the number of mismatches allowed. The resulting trie and index is written to a file with the input whtelist filename as the prefix, and `k{kmer_length}.pkl.gz` as the suffix.

```bash
scarecrow encode --barcodes BC1:v1:v1_whitelist.txt --pickle --index kmer --mismatches 2
```
