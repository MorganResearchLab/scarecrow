# scarecrow

<img style="float:right;width:200px;" src="./img/scarecrow.png" alt="scarecrow"/>

A library-agnostic toolkit for pre-processing combinatorial indexed single-cell RNA sequence data.

[Documentation](docs/root.md)

**scarecrow is undergoing substantial editing and may not behave as intended.**

### Issues

* Change pngs from harvest to output to barcode_positions basename rather than log basename
* FASTQ to SAM might be handy
* Input validation needs coding
* Jitter does not currently apply to UMI or insert sequence
  - if UMI on same read and downstream then may need position updating before extraction
* Docs need updating
  - revised SAM output
  - examples