# scarecrow

<img style="float:right;width:200px;" src="./img/scarecrow.png" alt="scarecrow"/>

**scarecrow is being actively developed, some kinks still being ironed out, and so it may not currently behave as expected.**

A library-agnostic toolkit for pre-processing combinatorial indexed single-cell RNA sequence data.

* [Installation](./docs/installation.md)
* [Workflow](./docs/workflow.md)
* [Troubleshooting](./docs/troubleshooting.md)
* Examples
    * [Parse Evercode WTv2](./docs/example_evercode.md)
    * [Scale Bio RNA](./docs/example_scale.md)
* Toolkit
    * [`encode`](./docs/toolkit_encode.md)
    * [`harvest`](./docs/toolkit_harvest.md)
    * [`json`](./docs/toolkit_json.md)
    * [`reap`](./docs/toolkit_reap.md)
    * [`recast`](./docs/toolkit_recast.md)
    * [`seed`](./docs/toolkit_seed.md)
    * [`sift`](./docs/toolkit_sift.md)
    * [`stats`](./docs/toolkit_stats.md)
    * [`weed`](./docs/toolkit_weed.md)

### To-do
* Ensure Aho-Corasick Trie matching method returns same results as set-based matching method
  - If pickle file is missing gz suffix then it needs adding automatically
  - Position counts returned are different vs set matched results
  - Unit-testing with reads generated to check edge cases
* Check/fix harvest image file indices (Scale returns 0,1,2,2,3,3; Parse returns 1,1,2,2)
* Jitter does not currently apply to UMI or insert sequence
  - if UMI on same read as barcodes, and downstream of any, then position *may* need adjusting dynamically before extracting sequence and qualities

### Wishlist
* Some tools would benefit from multi-threading (i.e. recast, stats)
* Weed should be refactored to work on FASTQ as well as SAM
