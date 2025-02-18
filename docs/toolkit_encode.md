<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[root](root.md)

## scarecrow encode
This tool generates an Aho-Corasick trie from a whitelist of barcodes for **exact** matching millions of short strings in parallel. Generating a trie for 6.9M 16-mer barcodes, pickling and compressing, takes around 3m30s. Note, the barcode label is stored in the trie, so the same trie should not be shared between barcodes even if based on the same whitelist.

```bash
scarecrow encode --barcodes BC1:v1:<v1_whitelist.txt> --out <BC1_v1_whitelist.trie.gz>    
```
