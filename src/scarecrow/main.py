#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scarecrow

A toolkit for pre-processing single-cell sequencing data.

"""
__author__ = "David Wragg"
__license__ = "GNU GPL v3.0"

from . import __version__
import sys
import argparse
import warnings
from scarecrow.seed import parser_seed, validate_seed_args
from scarecrow.reap import parser_reap, validate_reap_args
from scarecrow.harvest import parser_harvest, validate_harvest_args

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
        "seed": parser_seed(subparsers),
        "reap": parser_reap(subparsers),
        "harvest": parser_harvest(subparsers)
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
        "seed": validate_seed_args,
        "reap": validate_reap_args,
        "harvest": validate_harvest_args
    }
    COMMAND_TO_FUNCTION[sys.argv[1]](parser, args)


if __name__ == "__main__":
    main()    

