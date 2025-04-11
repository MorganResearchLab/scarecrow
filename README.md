# scarecrow

<img style="float:right;width:200px;" src="./img/scarecrow.png" alt="scarecrow"/>

A library-agnostic toolkit for pre-processing combinatorial indexed single-cell RNA sequence data.

[Documentation](docs/root.md)

**scarecrow is undergoing substantial editing and may not behave as intended.**

### Issues

* FASTQ output
  - corrected barcode quality scores are written, with N quality for negative space bases
  - this may cause an issue with some software as N is not a Phred quality score
* Input validation needs coding
* Jitter does not currently apply to UMI or insert sequence
  - if UMI on same read and downstream then may need position updating before extraction