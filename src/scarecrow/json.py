# -*- coding: utf-8 -*-
"""
#!/usr/bin/env python3
@author: David Wragg
"""

import json
import gzip
import logging
from pathlib import Path
from argparse import RawTextHelpFormatter
from scarecrow import __version__
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string
from scarecrow.recast import parse_tags_from_header


def parser_json(parser):
    subparser = parser.add_parser(
        "json",
        description="""
Generate a JSON file to accompany a scarecrow FASTQ file.

Example:

scarecrow json --sam cdna.sam
---
""",
        help="Converts a SAM file to a scarecrow interleaved FASTQ with accompanying JSON file.",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument(
        "-i",
        "--in",
        metavar="<file>",
        dest="infile",
        help=("FASTQ file to generate JSON"),
        type=str,
        required=True,
        default=[],
    )
    return subparser


def validate_json_args(parser, args) -> None:
    """
    Validate arguments
    """
    # Global logger setup
    logfile = "{}_{}.{}".format(
        "./scarecrow_json", generate_random_string(), "log"
    )
    logger = setup_logger(logfile)
    logger.info(f"scarecrow version {__version__}")
    logger.info(f"logfile: '{logfile}'")

    run_json(infile=args.infile)


@log_errors
def run_json(infile: str = None) -> None:
    """
    Main function to extract sequences and barcodes
    """
    logger = logging.getLogger("scarecrow")

    # Validate file exists
    if infile:
        if Path(infile).exists():
            json_file = infile.replace(".fastq", ".json")
            json_file = json_file[:-3] if json_file.endswith('.gz') else json_file
            logger.info(f"Generating JSON file '{json_file}' from '{infile}'")
            generate_json_from_fastq(fastq_file=infile, json_file=json_file)

        else:
            logger.info(f"'{infile}' does not exist")
    else:
        logger.info("No FASTQ file provided")

    logger.info("Finished!")



def generate_json_from_fastq(fastq_file: str, json_file: str) -> None:
    """
    Generate JSON file describing the FASTQ structure based on the first read's barcode and UMI information.

    Args:
        fastq_file: Path to the input FASTQ file
        json_file: Path to the output JSON file
    """
    logger = logging.getLogger("scarecrow")

    # Initialize variables
    barcode_lengths = []
    umi_length = 0
    whitelists = []

    # Open the FASTQ file (handles both gzipped and plain text)
    open_func = gzip.open if fastq_file.endswith('.gz') else open
    try:
        with open_func(fastq_file, 'rt') as fq:
            # Read first record
            header = fq.readline().strip()
            if not header:
                logger.error("Empty FASTQ file")
                return

            sequence = fq.readline().strip()
            _ = fq.readline()  # '+'
            _ = fq.readline()  # qualities

            # Parse tags from header
            tags = parse_tags_from_header(header)
            tags_dict = {tag: value for tag, value in tags}

            # Get barcode (CB tag) and UMI (UR tag)
            barcode = tags_dict.get("CB", "")
            umi = tags_dict.get("UR", "")

            # Calculate barcode lengths
            if barcode:
                barcode_parts = barcode.split(',')
                barcode_lengths = [len(part) for part in barcode_parts]
                # For whitelists, we'll use empty strings since we don't know the actual paths
                whitelists = [""] * len(barcode_parts)

            # Calculate UMI length
            if umi:
                umi_length = len(umi)

            # Verify lengths match the sequence
            total_expected_length = sum(barcode_lengths) + umi_length
            if total_expected_length > 0 and len(sequence) < total_expected_length:
                logger.warning(f"Sequence length ({len(sequence)}) is shorter than expected barcode+UMI length ({total_expected_length})")

    except Exception as e:
        logger.error(f"Error reading FASTQ file: {str(e)}")
        return

    # Generate JSON data
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
            "range": f"1:{current_position + 1}-{end_position}",
            "whitelist": whitelists[i] if i < len(whitelists) else ""
        })

        if kb_x is None:
            kb_x = f"0,{current_position},{end_position}"
            star_x = f"0_{current_position}_0_{end_position}"
        else:
            kb_x = f"{kb_x},0,{current_position},{end_position}"
            star_x = f"{star_x} 0_{current_position}_0_{end_position}"
        current_position = end_position

    # UMI information if present
    #star_umi = None
    if umi_length > 0:
        json_data["umi"].append({
            "range": f"1:{current_position + 1}-{current_position + umi_length}"
        })
        kb_x = f"{kb_x}:0,{current_position},{current_position + umi_length}"
        #star_umi = f"0_{current_position},0,{current_position + umi_length}"

    # Add kallisto-bustools command template
    if kb_x:
        json_data["kallisto-bustools"].append({
            "kb count": f"-i </path/to/transcriptome.idx> -g </path/to/transcripts_to_genes> -x {kb_x}:1,0,0 -w NONE --h5ad --inleaved -o <outdir> {fastq_file}"
        })

    # Write JSON file
    try:
        with open(json_file, "w") as f:
            json.dump(json_data, f, indent=4)
            f.write('\n')
        logger.info(f"Successfully generated JSON file: {json_file}")
    except Exception as e:
        logger.error(f"Error writing JSON file: {str(e)}")
