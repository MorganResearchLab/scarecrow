<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[root](root.md)

## scarecrow sift
The `sift` tool reads a SAM file (--sam) or a scarecrow interleaved FASTQ file (--fastq) and accompanying JSON file (--json), and filters out reads containing any invalid barcodes. An invalid barcode in a SAM file is a corrected barcode (CB) that contains an `N`. Reads without the CB tag are retained in the output. An invalid barcode in an interleaved FASTQ file is any barcode sequence range identified from the JSON file found to contain an `N`. The results are written to a new file with the `_sift` in the suffix (.e.g `_sift.sam` or `_sift.fastq`). 

```bash
scarecrow samstat --sam file.sam
```

