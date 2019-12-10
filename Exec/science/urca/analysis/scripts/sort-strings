#!/usr/bin/env python
""" 
Sort the input filenames and print them.

Donald E. Willcox
"""
import argparse
from FileToolkit import sort_input_filenames

parser = argparse.ArgumentParser()
parser.add_argument('infiles', type=str, nargs='+',
                    help='List of input files for sorting.')
args = parser.parse_args()


if __name__=='__main__':
    # Sort input files if needed
    input_files = sort_input_filenames(args.infiles)

    # Print sorted files
    for f in input_files:
        print(f)
