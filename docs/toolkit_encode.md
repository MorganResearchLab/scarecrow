<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[root](../README.md)

## scarecrow encode
This tool generates an Aho-Corasick trie and k-mer index from a whitelist of barcodes. Generating a trie for 6.9M 16-mer barcodes, pickling and compressing, takes around 3m30s. Note, the barcode label is stored in the trie, so the same trie should not be shared between barcodes even if based on the same whitelist. The trie is used for exact matching, the k-mer index is used for approximate matching to reduce the search space. The trie and k-mer index is written to a file with the input whtelist filename as the prefix, and `k{kmer_length}.pkl.gz` as the suffix.

The `kmer_length` determines the size of the k-mer when generating the index and this has a functional impact on approximate matching when considering mismatches with `scarecrow reap`. For example, if the barcode length is 8 and a k-mer length of 4 is used, then barcode index ACGTTGCA will have the k-mer values ACGT, CGTT, GTTG, TTGC, TGCA. If we then consider 2 mismatches, i.e. ACaTTaCA, then these (lowercase) mismatches impact all k-mers of the barcode index and no match will be found. At most, to ensure at least one error-free k-mer, the maximum k-mer length should be set to
$$
m k < n - k + 1 \;\;\Longrightarrow\;\; k \le \frac{n}{m+1}
$$
where $$ n $$ is the barcode length and $$ m $$ is the number of mismatches allowed.
$$
k_{\max} = \left\lfloor \tfrac{n}{m+1} \right\rfloor
$$

```bash
scarecrow encode --barcodes BC1:v1:v1_whitelist.txt --pickle --kmer_length 8
```
