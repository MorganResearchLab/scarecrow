# scarecrow

<img style="float:right;width:200px;" src="./img/scarecrow.png" alt="scarecrow"/>

**scarecrow is being actively developed, some kinks still being ironed out, and so it may not currently behave as expected.**

A library-agnostic toolkit for pre-processing combinatorial indexed single-cell RNA sequence data.

[Documentation](docs/root.md)

### Issues

* Docs need updating (esp examples)
* Possible issue with harvest returning 1-based instead of 0-based file index in some instances
  - WTv2 test data gives barcode_positions_set.csv file_index 2 (reap returns error, but works with file_index 1)
* Jitter does not currently apply to UMI or insert sequence
  - if UMI on same read as barcodes, and downstream of any, then position may need adjusting dynamically before extracting sequence and qualities

### Wishlist

* Some tools would benefit from multi-threading (i.e. recast, stats)
