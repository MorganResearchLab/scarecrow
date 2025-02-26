<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[root](root.md)

## scarecrow encode
This tool generates an Aho-Corasick trie and k-mer index from a whitelist of barcodes. Generating a trie for 6.9M 16-mer barcodes, pickling and compressing, takes around 3m30s. Note, the barcode label is stored in the trie, so the same trie should not be shared between barcodes even if based on the same whitelist. The trie is used for exact matching, the k-mer index is used for approximate matching to reduce the search space. 

The `kmer_length` determines the size of the k-mer when generating the index and this has a functional impact on approximate matching when considering mismatches with the `reap` tool. For example, an if the barcode length is 8 and a k-mer length of 4 is used, then the barcode of ACGTTGCA will have the values ACGT, CGTT, GTTG, TTGC, TGCA in the index. If we then consdider 2 mismatches, i.e. ACaTTaCA, then these mismatches impact all k-mers for the barcode in the index and no match will be found. A smaller k-mer length, i.e. 2, would be required in this instance.

```bash
scarecrow encode --barcodes BC1:v1:v1_whitelist.txt --trie --kmer_length 4
```
