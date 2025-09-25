# -*- coding: utf-8 -*-
"""
#!/usr/bin/env python3
@author: David Wragg
"""

import re
import gzip
import json
import pysam
import logging
from pathlib import Path
from argparse import RawTextHelpFormatter
from scarecrow import __version__
from scarecrow.logger import log_errors, setup_logger
from scarecrow.tools import generate_random_string
from dataclasses import fields


def parser_recast(parser):
    subparser = parser.add_parser(
        "recast",
        description="""
Recasts either a scarecrow SAM file to an interleaved FASTQ with accompanying JSON file, or a scarecrow FASTQ file into a SAM file.

Example:

scarecrow recast --in cdna.sam
---
""",
        help="Recasts either a SAM file to a scarecrow interleaved FASTQ with accompanying JSON file, or a FASTQ file into a SAM file.",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument(
        "-i",
        "--in",
        metavar="<file>",
        dest="infile",
        help=("File to recast"),
        type=str,
        required=True,
        default=[],
    )
    return subparser


def validate_recast_args(parser, args) -> None:
    """
    Validate arguments
    """
    # Global logger setup
    logfile = "{}_{}.{}".format(
        "./scarecrow_recast", generate_random_string(), "log"
    )
    logger = setup_logger(logfile)
    logger.info(f"scarecrow version {__version__}")
    logger.info(f"logfile: '{logfile}'")

    run_recast(infile=args.infile, args_string=" ".join(f"--{k} {v}" for k, v in vars(args).items() if v is not None))


@log_errors
def run_recast(infile: str = None, args_string: str = None) -> None:
    """
    Main function to extract sequences and barcodes
    """
    logger = logging.getLogger("scarecrow")

    # Determine input file type
    input_path = Path(infile)
    if not input_path.exists():
        logger.error(f"Input file '{infile}' does not exist")
        return

    if infile.lower().endswith(('.sam', '.bam')):
        # SAM to FASTQ conversion
        fastq_file = input_path.with_suffix('.fastq').as_posix()
        json_file = input_path.with_suffix('.json').as_posix()
        logger.info(f"Recasting '{infile}' to interleaved FASTQ '{fastq_file}'")
        logger.info(f"Will generate JSON file: '{json_file}'")
        run_sam2fastq(sam_file=infile, fastq_file=fastq_file, json_file=json_file)

    elif infile.lower().endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
        # FASTQ to SAM conversion
        base_path = input_path
        if infile.lower().endswith('.gz'):
            base_path = input_path.with_suffix('')  # Remove .gz
        output_sam = base_path.with_suffix('.sam').as_posix()
        logger.info(f"Converting FASTQ '{infile}' to SAM '{output_sam}'")
        run_fastq2sam(fastq_file=infile, output_sam=output_sam, args_string=args_string)

    else:
        logger.error(f"Unrecognized file format for '{infile}'. Must be .sam or .fastq")
        return

    logger.info("Finished!")

@log_errors
def run_sam2fastq(sam_file: str = None, fastq_file: str = None, json_file: str = None) -> None:
    """
    Convert SAM to interleaved FASTQ with JSON metadata
    """
    logger = logging.getLogger("scarecrow")

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
                # Split barcode
                barcode_parts = barcode.split('_')
                barcode_lengths = [len(part) for part in barcode_parts]

            if umi and umi_length is None:
                umi_length = len(umi)

            # Construct comprehensive header with all tags
            header_parts = []
            for tag, value in read.tags:
                header_parts.append(f"{tag}={value}")
            header = f"@{read.query_name} " + " ".join(header_parts)

            # R1 (header contains all tags)
            r1_header = f"{header} /1"
            r1_seq = barcode.replace("_", "")
            r1_qual = tags.get("XQ", "F" * len(r1_seq))
            if umi:
                r1_seq += umi
                r1_qual += tags.get("UY", "F" * umi_length)

            # R2 (sequence from SAM)
            r2_header = f"{header} /2"
            r2_seq = read.query_sequence
            r2_qual = "".join(chr(q + 33) for q in read.query_qualities)

            # Write interleaved FASTQ
            fq.write(f"{r1_header}\n{r1_seq}\n+\n{r1_qual.replace("_", "")}\n")
            fq.write(f"{r2_header}\n{r2_seq}\n+\n{r2_qual}\n")

        # Generate JSON file
        generate_json(
            barcode_lengths=barcode_lengths,
            umi_length=umi_length,
            json_file=json_file,
            fastq_file=fastq_file
        )


@log_errors
def run_fastq2sam(fastq_file: str = None, output_sam: str = None, args_string: str = None) -> None:
    """
    Convert interleaved FASTQ (with tags in header) back to SAM format
    """
    logger = logging.getLogger("scarecrow")

    # Prepare SAM header (minimal by default)
    escaped_cmd = args_string.replace("\t", " ").replace("\n", " ")
    header = {'HD': {'VN': '1.6'},
              'PG': [{'ID': 'reap',
                     'PN': 'scarecrow',
                     'VN': __version__,
                     'CL': escaped_cmd}] }

    open_func = gzip.open if fastq_file.endswith('.gz') else open
    with pysam.AlignmentFile(output_sam, "wh", header=header) as sam_out:
        with open_func(fastq_file, 'rt') as fq:
            while True:
                # Read R1 (header contains all tags)
                r1_header = fq.readline().strip()
                if not r1_header:  # End of file
                    break

                # Skip R1 sequence and qualities (not needed)
                _ = fq.readline()  # Sequence
                _ = fq.readline()  # '+'
                _ = fq.readline()  # Qualities

                # Read R2 (actual sequence to include in SAM)
                r2_header = fq.readline().strip()
                r2_seq = fq.readline().strip()
                _ = fq.readline()  # '+'
                r2_qual = fq.readline().strip()

                # Extract read name (before first space)
                read_name = r1_header.split(' ')[0][1:]  # Remove @

                # Parse all tags from header
                tags = parse_tags_from_header(r1_header)

                # Create SAM record
                record = pysam.AlignedSegment()
                record.query_name = read_name
                record.query_sequence = r2_seq
                record.query_qualities = pysam.qualitystring_to_array(r2_qual)
                record.flag = 4  # unmapped
                record.tags = tags

                sam_out.write(record)

def parse_tags_from_header(header: str) -> list:
    """
    Parse SAM tags from FASTQ header string.

    Header format example:
    @SRR28867558.10001 CR=... CY=... CB=... XQ=... XP=... XM=... UR=... UY=... /1

    Returns:
        List of (tag, value) pairs for SAM record
    """
    #logger = logging.getLogger("scarecrow")
    tags = []
    fields = header.split()
    #logger.info(f"{header}")
    allowed_keys = {"CB", "CR", "CY", "XQ", "XP", "XM", "UY", "UR"}
    for field in fields:
        if '=' in field:
            key, _, value = field.partition('=')
            if key in allowed_keys:
                # Remove trailing /<number> if present
                value = re.sub(r'/\d+$', '', value)
                tags.append((key,value))
    #logger.info(f"{tags}")

    return tags


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
            "range": f"1:{current_position + 1}-{end_position}",
            "whitelist": ""  # Empty since we don't have whitelist info from SAM
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
    if umi_length is not None:
        json_data["umi"].append({
            "range": f"1:{current_position + 1}-{current_position + umi_length}"
        })
        kb_x = f"{kb_x}:0,{current_position},{current_position + umi_length}"
        #star_umi = f"0_{current_position},0,{current_position + umi_length}"

    # Add kallisto-bustools command template
    json_data["kallisto-bustools"].append({
        "kb count": f"-i </path/to/transcriptome.idx> -g </path/to/transcripts_to_genes> -x {kb_x}:1,0,0 -w NONE --h5ad --inleaved -o <outdir> {fastq_file}"
    })

    # Write JSON file
    with open(json_file, "w") as f:
        json.dump(json_data, f, indent=4)
        f.write('\n')
