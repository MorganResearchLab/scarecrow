# -*- coding: utf-8 -*-
"""
#!/usr/bin/env python3
@author: David Wragg
"""

import json
import pysam
import logging
from pathlib import Path
from argparse import RawTextHelpFormatter
from scarecrow import __version__
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string


def parser_sam2fastq(parser):
    subparser = parser.add_parser(
        "sam2fastq",
        description="""
Convert a SAM file to a scarecrow interleaved FASTQ with accompanying JSON file.

Example:

scarecrow sam2fastq --sam cdna.sam
---
""",
        help="Converts a SAM file to a scarecrow interleaved FASTQ with accompanying JSON file.",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument(
        "-s",
        "--sam",
        metavar="<file>",
        help=("SAM file to convert to interleaved FASTQ"),
        type=str,
        required=True,
        default=[],
    )
    return subparser


def validate_sam2fastq_args(parser, args) -> None:
    """
    Validate arguments
    """
    # Global logger setup
    logfile = "{}_{}.{}".format(
        "./scarecrow_sam2fastq", generate_random_string(), "log"
    )
    logger = setup_logger(logfile)
    logger.info(f"scarecrow version {__version__}")
    logger.info(f"logfile: '{logfile}'")

    run_sam2fastq(sam_file=args.sam)


@log_errors
def run_sam2fastq(sam_file: str = None) -> None:
    """
    Main function to extract sequences and barcodes
    """
    logger = logging.getLogger("scarecrow")

    # Validate file exists
    if sam_file:
        if Path(sam_file).exists():
            fastq_file = sam_file.replace(".sam", ".fastq")
            json_file = sam_file.replace(".sam", ".json")
            logger.info(f"Converting '{sam_file}' to interleaved FASTQ '{fastq_file}'")
            logger.info(f"Will generate JSON file: '{json_file}'")

            # Track barcode and UMI lengths for JSON generation
            barcode_lengths = []
            umi_length = None

            with (
                pysam.AlignmentFile(sam_file, "rb", check_sq=False) as sam,
                open(fastq_file, "w") as fq,
            ):
                # Force reading even if header is missing
                for read in sam.fetch(until_eof=True):

                    # Extract tags
                    tags = {k: str(v) for k, v in read.tags}
                    
                    # Get barcode (CB tag) and UMI (UR tag)
                    barcode = tags.get("CB", "")
                    umi = tags.get("UR", "")
                    
                    # Track lengths for JSON generation
                    if barcode and not barcode_lengths:
                        # Split barcode by underscores
                        barcode_parts = barcode.split('_')
                        barcode_lengths = [len(part)-1 for part in barcode_parts]
                    
                    if umi and umi_length is None:
                        umi_length = len(umi)

                    # R1 (barcode + UMI)
                    r1_header = f"@{read.query_name}/1"
                    r1_seq = barcode.replace("_", "")
                    r1_qual = tags.get("CY", "F" * len(r1_seq)).replace("_", "")  # Default to high quality if no quality scores
                    if umi:
                        r1_seq += umi                    
                        r1_qual += tags.get("UY", "F" * umi_length).replace("_", "")  # Default to high quality if no quality scores
                    
                    # R2 (sequence from SAM)
                    r2_header = f"@{read.query_name}/2"
                    r2_seq = read.query_sequence
                    r2_qual = "".join(chr(q + 33) for q in read.query_qualities)
                    
                    # Write interleaved FASTQ
                    fq.write(f"{r1_header}\n{r1_seq}\n+\n{r1_qual}\n")
                    fq.write(f"{r2_header}\n{r2_seq}\n+\n{r2_qual}\n")

                # Generate JSON file
                generate_json(
                    barcode_lengths = barcode_lengths,
                    umi_length = umi_length,
                    json_file = json_file,
                    fastq_file = fastq_file
                )


        else:
            logger.info(f"'{sam_file}' does not exist")
    else:
        logger.info("No SAM file provided")

    logger.info("Finished!")

def generate_json(barcode_lengths: list, umi_length: int, json_file: str, fastq_file: str) -> None:
    """
    Generate JSON file describing the FASTQ structure based on barcode and UMI lengths
    """
    json_data = {
        "description": "scarecrow",
        "barcodes": [],
        "umi": [],
        "kallisto-bustools": []
    }
    
    # Barcode information
    current_position = 0
    kb_x = None
    star_x = None
    
    for i, length in enumerate(barcode_lengths):
        end_position = current_position + length
        json_data["barcodes"].append({
            "range": f"1:{current_position}-{end_position}",
            "whitelist": ""  # Empty since we don't have whitelist info from SAM
        })
        
        if kb_x is None:
            kb_x = f"0,{current_position},{end_position}"
            star_x = f"0_{current_position}_0_{end_position}"
        else:
            kb_x = f"{kb_x},0,{current_position},{end_position}"
            star_x = f"{star_x} 0_{current_position}_0_{end_position}"
        current_position = end_position + 1
    
    # UMI information if present
    star_umi = None
    if umi_length is not None:
        json_data["umi"].append({
            "range": f"1:{current_position}-{current_position + umi_length - 1}"
        })
        kb_x = f"{kb_x}:0,{current_position},{current_position + umi_length - 1}"
        star_umi = f"0_{current_position},0,{current_position + umi_length - 1}"
    
    # Add kallisto-bustools command template
    json_data["kallisto-bustools"].append({
        "kb count": f"-i </path/to/transcriptome.idx> -g </path/to/transcripts_to_genes> -x {kb_x}:1,0,0 -w NONE --h5ad --inleaved -o <outdir> {fastq_file}"
    })

    # Write JSON file
    with open(json_file, "w") as f:
        json.dump(json_data, f, indent=4)

