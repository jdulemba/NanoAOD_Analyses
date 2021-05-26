from coffea.util import load, save
import re
from pdb import set_trace
import coffea.processor as processor
import coffea.hist

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('dataset_to_sub', type=str, help='Name of dataset to be substituted from new_file to replace in orig_file.')
parser.add_argument('output_fname', type=str, help='Name of output file with file extension.')
parser.add_argument('orig_file', type=str, help="Input files separated by ':'")
parser.add_argument('new_file', type=str, help="Input files separated by ':'")
args = parser.parse_args()

orig_dict = load(args.orig_file)
new_dict = load(args.new_file)
output_dict = orig_dict.copy()

    # get mask of all datasets that will remain unchanged
non_dataset_mask = re.compile('(?!%s*)' % args.dataset_to_sub)
dataset_mask = re.compile('(%s*)' % args.dataset_to_sub)

for hname in new_dict.keys():
    if hname not in orig_dict.keys():
        print('%s is in new_dict but not orig_dict, skipping')
        continue
    if 'cutflow' in hname: continue
    print(hname)
    new_obj = new_dict[hname]
    orig_obj = orig_dict[hname]

    if isinstance(new_obj, coffea.hist.hist_tools.Hist) and isinstance(orig_obj, coffea.hist.hist_tools.Hist):
        out_histo = output_dict[hname].copy() # make copy of output histogram
        out_histo.clear() # clear output histogram
        out_histo.add(new_obj[dataset_mask]) # add substituted histogram to output
        out_histo.add(orig_obj[non_dataset_mask]) # add other histograms from original histogram to output
        output_dict[hname] = out_histo # substitue new histogram into output_dict

    elif isinstance(new_obj, processor.accumulator.defaultdict_accumulator) and isinstance(orig_obj, processor.accumulator.defaultdict_accumulator):
        output_dict[hname] = new_obj

    elif isinstance(new_obj, processor.accumulator.value_accumulator) and isinstance(orig_obj, processor.accumulator.value_accumulator):
        output_dict[hname] = new_obj

    else:
        raise ValueError("type(new_obj) %s and type(orig_obj) %s are not currently supported or are different." % (type(new_obj), type(orig_obj)))

#set_trace()
outname = args.output_fname if args.output_fname.endswith('.coffea') else '%s.coffea' % args.output_fname
output_acc = processor.dict_accumulator(output_dict)
save(output_acc, outname)
print('%s written' % outname)
