import gzip
from Bio.Seq import Seq
from collections import defaultdict
from argparse import RawTextHelpFormatter

def parser_barcode_check(parser):
    subparser = parser.add_parser(
        "barcode_check",
        description="""
Check barcode whitelist against a fastq file and return counts of matched barcode sequences and their start position

Examples:
scarecrow check_barcodes barcodes.txt R1.fastq.gz barcode_counts.txt
---
""",
        help="Check barcodes",
        formatter_class=RawTextHelpFormatter,
    )
    subparser.add_argument("barcodes", help="Text file list of barcodes (one per line)")
    subparser.add_argument("fastq", help="A single FASTQ files")
    subparser.add_argument("counts", help="Output file for barcode counts")
    return subparser

def validate_barcode_check_args(parser, args):
    run_barcode_check(barcode_file = args.barcodes, fastq_file = args.fastq, counts_file = args.counts)

def run_barcode_check(barcode_file, fastq_file, counts_file):
    # Load barcodes
    barcodes = load_barcodes(barcode_file)
    # Count barcodes
    barcode_data = count_barcodes(fastq_file, barcodes)
    # Save results
    save_results(barcode_data, counts_file)
    print(f"Results saved to {counts_file}")


def load_barcodes(barcode_file):
    """
    Load barcodes from a file into a set.
    """
    with open(barcode_file, 'r') as f:
        barcodes = {line.strip() for line in f}
    return barcodes

def count_barcodes(fastq_file, barcodes):
    """
    Count occurrences of barcodes and their reverse complements in a FASTQ file.
    Track positions of matches in the sequences.
    """
    barcode_data = defaultdict(lambda: {"direct": 0, "reverse": 0, "positions": []})
    reverse_complements = {barcode: str(Seq(barcode).reverse_complement()) for barcode in barcodes}
    
    with gzip.open(fastq_file, 'rt') as f:
        for i, line in enumerate(f):
            # Process only sequence lines (2nd line in every 4-line FASTQ entry)
            if i % 4 == 1:
                sequence = line.strip()
                for barcode in barcodes:
                    # Check forward barcode
                    pos = sequence.find(barcode)
                    while pos != -1:
                        barcode_data[barcode]["direct"] += 1
                        barcode_data[barcode]["positions"].append(f"F:{pos}")
                        pos = sequence.find(barcode, pos + 1)

                    # Check reverse complement
                    rev_barcode = reverse_complements[barcode]
                    pos = sequence.find(rev_barcode)
                    while pos != -1:
                        barcode_data[barcode]["reverse"] += 1
                        barcode_data[barcode]["positions"].append(f"R:{pos}")
                        pos = sequence.find(rev_barcode, pos + 1)

    return barcode_data

def save_results(barcode_data, output_file):
    """
    Save barcode counts, direction annotations, and positions to a text file.
    """
    with open(output_file, 'w') as f:
        f.write("Barcode\tCount\tDirection\tPositions\n")
        for barcode, data in barcode_data.items():
            total_count = data["direct"] + data["reverse"]
            if data["direct"] > data["reverse"]:
                direction = "forward"
            elif data["reverse"] > data["direct"]:
                direction = "reverse complement"
            else:
                direction = "both"  # In case counts are equal
            positions = ",".join(data["positions"]) if data["positions"] else "None"
            f.write(f"{barcode}\t{total_count}\t{direction}\t{positions}\n")


