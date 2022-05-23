#!/usr/bin/env python

import Utilities.plot_tools as plt_tools
from pdb import set_trace

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('output_fname', type=str, help='Name of output file with file extension.')
parser.add_argument('input_files', type=str, help="Input files separated by ':'")
parser.add_argument("--add_dict", action="store_true", help="Inputs are dictionaries.")
args = parser.parse_args()

input_files = args.input_files.split(':')


if args.add_dict:
    import coffea.processor.accumulator as proc_acc
    from coffea.util import load
    from tqdm import tqdm

    if len(input_files) < 2:
        raise ValueError("Input file list must have at least 2 inputs")

    first_acc = load(input_files[0]).copy()
    second_acc = load(input_files[1]).copy()

    tot_dists = 0    
    print("0:", len(first_acc["3Jets"]["Muon"].keys()))
    tot_dists += len(first_acc["3Jets"]["Muon"].keys())
    print("1:", len(second_acc["3Jets"]["Muon"].keys()))
    tot_dists += len(second_acc["3Jets"]["Muon"].keys())

    #set_trace()
    output_acc = proc_acc.add(first_acc, second_acc)
    for idx in tqdm(range(len(input_files))):
        if idx < 2:
            continue
        else:
            try:
                tmp_acc = load(input_files[idx]).copy()
                print(f"{idx}:", len(tmp_acc["3Jets"]["Muon"].keys()))
                tot_dists += len(tmp_acc["3Jets"]["Muon"].keys())
                output_acc = proc_acc.add(output_acc, load(input_files[idx]).copy() )
            except:
                raise ValueError(f"File {input_files[idx]} (number {idx}) could not be added")
    print("Final output:", len(output_acc["3Jets"]["Muon"].keys()))
    print("Final expected:", tot_dists)
    #set_trace()
    plt_tools.save_accumulator(output_acc, args.output_fname)

else:
    output_acc = plt_tools.add_coffea_files(input_files)
    plt_tools.save_accumulator(output_acc, args.output_fname)

