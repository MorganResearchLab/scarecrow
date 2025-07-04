# scarecrow

<img style="float:right;width:200px;" src="./img/scarecrow.png" alt="scarecrow"/>

**scarecrow is being actively developed, some kinks still being ironed out, and so it may not currently behave as expected.**

A library-agnostic toolkit for pre-processing combinatorial indexed single-cell RNA sequence data.

[Documentation](docs/root.md)

### Issues

* Harvest image files from Scale have indices 0,1,2,2,3,3 whilst Parse correctly show 1,1,2,2
* Update docs
* Unit tests
* Jitter does not currently apply to UMI or insert sequence
  - if UMI on same read as barcodes, and downstream of any, then position may need adjusting dynamically before extracting sequence and qualities

### Wishlist

* Some tools would benefit from multi-threading (i.e. recast, stats)
* Weed should be refactored to work on FASTQ as well as SAM
