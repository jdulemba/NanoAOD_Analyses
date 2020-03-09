#!/usr/bin/env python

from coffea.util import load, save
from pdb import set_trace
import collections
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('output_fname', type=str, help='Name of output file with file extension.')
parser.add_argument('input_files', type=str, help="Input files separated by ':'")

args = parser.parse_args()

input_files = args.input_files.split(':')
input_accs = [load(fname) for fname in input_files]

output_acc = collections.Counter()
for acc in input_accs:
    output_acc.update(acc)

save(output_acc, args.output_fname)
print('%s written' % args.output_fname)
