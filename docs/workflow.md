<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[Back to root](root.md)

# Workflow

scarecrow is a scRNA-seq pre-processing tool developed to increase the usable data generated from an experiment by idenifying jittered barcodes, correcting mismatches, and outputting the data in a format used by popular third-party tools. 

## Pseudoaligment

If the aim is to perform cell-gene quantification via pseudoalignment using kallisto-bustools, then an interleaved FASTQ file should be generated with `scarecrow reap`. The interleaved FASTQ contains corrected barcodes followed by the UMI in read 1 and the extracted sequence in read 2, and is supplemented by a JSON file that contains the `kb count` parameters required to process the FASTQ file.

## Alignment

If the aim is to perform cell-gene quantification via standard alignment, for example the UMI-tools workflow, then a SAM file should be generated with `scarecrow reap`. The SAM file can be aligned with STAR, annoated with featureCounts, and processed with UMI-tools (taking advantage of the read tags in the SAM file) to generate a counts matrix.