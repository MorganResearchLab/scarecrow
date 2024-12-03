#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scarecrow

A toolkit to parse seqspec.yaml files for downstream analysis of single-cell sequencing data.

"""
__author__ = "David Wragg"
__license__ = "GNU GPL v3.0"

from . import __version__
import sys
import argparse
import warnings
from .scarecrow_extract import parser_extract, validate_extract_args
from .scarecrow_barcodes import parser_barcodes, validate_barcodes_args

def main():
    warnings.simplefilter("default", DeprecationWarning)
    
    # setup parsers
    parser = argparse.ArgumentParser(
        description=f"""
\033[38;5;202m
 _,  _,_   ,_   _, _,,_   _, ,  , 
(_, / '|\\  |_) /_,/  |_) / \\,| ,| 
 _)'\\_ |-\\'| \\'\\_'\\_'| \\'\\_/ |/\\| 
'     `'  `'  `  `  `'  `'   '  ` 
{__version__}
\033[32m
A toolkit to parse seqspec.yaml files for downstream analysis of single-cell sequencing data.
\033[0m
GitHub: https://github.com/MorganResearchLab/scarecrow
Documentation: https://www.morganlab.co.uk/software/scarecrow

""",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    subparsers = parser.add_subparsers(
        dest="command",
        metavar="<CMD>",
    )

    # Setup the arguments for all subcommands
    command_to_parser = {
        "extract": parser_extract(subparsers),
        "barcodes": parser_barcodes(subparsers)
    }
    
    # Show help when no arguments are given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    if len(sys.argv) == 2:
        if sys.argv[1] in command_to_parser:
            command_to_parser[sys.argv[1]].print_help(sys.stderr)
        elif sys.argv[1] == "--version":
            print(f"scarecrow {__version__}")
        else:
            parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()    

    # Setup validator and runner for all subcommands (validate and run if valid)
    COMMAND_TO_FUNCTION = {
        "extract": validate_extract_args,
        "barcodes": validate_barcodes_args
    }
    COMMAND_TO_FUNCTION[sys.argv[1]](parser, args)


if __name__ == "__main__":
    main()    

