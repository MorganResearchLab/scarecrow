# scarecrow

<img style="float:right;width:200px;" src="./img/scarecrow.png" alt="scarecrow"/>

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
    * [`inspect`](./docs/toolkit_inspect.md)
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
* FASTQ output sequence length is taken from one read but written in header for both
  - if R1 is len 28 and R2 is len 90, then in interleaved FASTQ read header output by scarecrow, both reads report the same length (e.g. 90)
  - need to store each read's sequence length and output the correct value
* Some tools would benefit from multi-threading (i.e. recast, stats)
* Would be helpful if `weed` could be refactored to work on FASTQ input format rather than just SAM
