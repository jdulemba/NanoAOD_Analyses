from coffea.util import load, save
import re
from pdb import set_trace
import coffea.processor as processor

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('output_fname', type=str, help='Name of output file with file extension.')
parser.add_argument('ttJets_PS_file', type=str, help="Input files separated by ':'")
parser.add_argument('ttJets_file', type=str, help="Input files separated by ':'")

args = parser.parse_args()

ttPS_dict = load(args.ttJets_PS_file)
tt_dict = load(args.ttJets_file)
output_dict = tt_dict.copy()

nonTTPS_mask = re.compile('(?!ttJets_PS*)')
sys_mask = re.compile('(?!nosys)')

for hname in ttPS_dict.keys():
    if 'cutflow' in hname: continue
    histo = output_dict[hname]
    #set_trace()
    ps_histo = ttPS_dict[hname]

    #    ## get ttJets_PS hists for nosys
    #tt_ps_nosys_dict = ps_histo['ttJets_PS*', 'nosys', :, :,:, :]
        ## get ttJets_PS hists for systematic variations
    tt_ps_sys_dict = ps_histo['ttJets_PS*', sys_mask, :, :, :, :]
        ## get all non ttJets_PS hists
    nonTTPS_dict = ps_histo[nonTTPS_mask, :, :, :, :, :]

        # add hists to ttJets hist
    histo.add(nonTTPS_dict)
    histo.add(tt_ps_sys_dict)


outname = args.output_fname if args.output_fname.endswith('.coffea') else '%s.coffea' % args.output_fname
output_acc = processor.dict_accumulator(output_dict)
save(output_acc, outname)
print('%s written' % outname)
