# scarecrow

<img style="float:right;width:200px;" src="./img/scarecrow.png" alt="scarecrow"/>

**scarecrow is subject to ongoing editing and may not behave as expected.**

A library-agnostic toolkit for pre-processing combinatorial indexed single-cell RNA sequence data.

[Documentation](docs/root.md)

### Issues

* Input validation needs coding (done: seed)
* Docs need updating (esp examples)
* Possible issue with harvest returning 1-based instead of 0-based file index in some instances
* Jitter does not currently apply to UMI or insert sequence
  - if UMI on same read as barcodes, and downstream of any, then position may need adjusting dynamically before extracting sequence and qualities

### Wishlist

* Some tools would benefit from multi-threading (i.e. recast, stats)