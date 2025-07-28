<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[Back to root](root.md)

# Workflow

scarecrow is a scRNA-seq pre-processing tool developed to increase the usable data generated from an experiment by idenifying jittered barcodes, correcting mismatches, and outputting the data in a format used by popular third-party tools. Below is a schematic diagram illustrating the scarecrow workflow. This begins with raw FASTQ files and barcode whitelists being passed to seed. The resulting barcode profiles are passed to harvest which outputs a file of expected barcode positions. These positions, together with the raw FASTQ files and barcode whitelists are passed to reap, which returns either a SAM or FASTQ file containing the target sequence and barcode information. Optional steps include extracting a barcode from the FASTQ read header with weed, if the reads were demultiplexed, and removing any reads containing unmatched barcodes with sift.

<img style="float:center;width:400px;" src="../img/scarecrow_flow.svg" alt="scarecrow flowchart"/>

We recommend adapter trimming after running scarecrow to ensure that read lengths are consistent when pre-processing the data, otherwise there is the potential for trimming to result in offset barcodes that might not subsequently be matched.


## Pseudoaligment

If the aim is to perform cell-gene quantification via pseudoalignment using kallisto-bustools, then an interleaved FASTQ file should be generated with `scarecrow reap`. The interleaved FASTQ contains corrected barcodes followed by the UMI in read 1 and the extracted sequence in read 2, and is supplemented by a JSON file that contains the `kb count` parameters required to process the FASTQ file. If it is necessary to add a barcode from the FASTQ header, then a SAM file should be generated with `scarecrow reap` and annotated with `scarecrow weed`. The SAM file can be converted to an interleaved FASTQ file for use with kallisto by `scarecrow recast`.

## Alignment

If the aim is to perform cell-gene quantification via standard alignment, for example the UMI-tools workflow, then a SAM file should be generated with `scarecrow reap`. The SAM file can be aligned with STAR, annoated with featureCounts, and processed with UMI-tools (taking advantage of the read tags in the SAM file) to generate a counts matrix.
