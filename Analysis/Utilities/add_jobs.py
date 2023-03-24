#!/usr/bin/env python

import Utilities.plot_tools as plt_tools
from pdb import set_trace

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("output_fname", type=str, help="Name of output file with file extension.")
parser.add_argument("input_files", type=str, help="Input files separated by ':'")
parser.add_argument("--add_dict", action="store_true", help="Inputs are dictionaries.")
args = parser.parse_args()

input_files = args.input_files.split(":")

if args.add_dict:
    output_acc = plt_tools.add_dict(input_files)
    plt_tools.save_accumulator(output_acc, args.output_fname)
else:
    output_acc = plt_tools.add_coffea_files(input_files)
    plt_tools.save_accumulator(output_acc, args.output_fname)

