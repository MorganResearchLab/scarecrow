# scarecrow

<img style="float:right;width:200px;" src="./img/scarecrow.png" alt="scarecrow"/>

A library-agnostic toolkit for pre-processing combinatorial indexed single-cell RNA sequence data.

[Documentation](docs/root.md)

**scarecrow is undergoing substantial editing and may not behave as intended.**

### Issues

* Possible issue with harvest returning 1-based instead of 0-based file index in some instances
* Check stats generated from reshaped SAM and FASTQ are consistent
* SAM (after weed) to FASTQ - need to check FASTQ and JSON output to ensure weeded barcode is present and correct
* Input validation needs coding  
* Jitter does not currently apply to UMI or insert sequence
  - if UMI on same read and downstream then may need position updating before extraction
* Docs need updating
  - revised SAM output
  - examples