# scarecrow

<img style="float:right;width:200px;" src="./img/scarecrow.png" alt="scarecrow"/>

**scarecrow is being actively developed, some kinks still being ironed out, and so it may not currently behave as expected.**

A library-agnostic toolkit for pre-processing combinatorial indexed single-cell RNA sequence data.

* [Installation](./docs/installation.md)
* [Introduction](./docs/workflow.md)
* [Troubleshooting](./docs/troubleshooting.md)
* Examples
    * [Parse Evercode WTv2](./docs/example_evercode.md)
    * [Scale QuantumnScale RNA](./docs/example_scale.md)
    * [Scarecrow paper code](https://github.com/MorganResearchLab/scarecrow_paper)
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
* Jitter does not currently apply to UMI or insert sequence
  - if UMI on same read as barcodes, and downstream of any, then position *may* need adjusting dynamically before extracting sequence and qualities

### Wishlist
* Some tools would benefit from multi-threading (i.e. recast, stats)
* Would be helpful if `weed` could be refactored to work on FASTQ input format rather than just SAM
