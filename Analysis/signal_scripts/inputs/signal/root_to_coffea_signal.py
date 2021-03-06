from coffea.lookup_tools.root_converters import convert_histo_root_file
from coffea.lookup_tools.dense_lookup import dense_lookup
from pdb import set_trace
from coffea.util import save
import os

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('input', help='Name of input root file')
parser.add_argument('--outname', help='Name of output coffea file')
args = parser.parse_args()

rfile = args.input if args.input.endswith('.root') else '%s.root' % args.input
if not os.path.isfile(rfile): raise ValueError('%s not found' % rfile)

signals = convert_histo_root_file(rfile)

output_dict = {}

topologies = sorted(list(set([key[0].split('/')[0] for key in signals.keys()])))
cvarname = 'topCosth_vs_mtt' # hardcoded

#varnames = list(set([key[0].split('/')[1] for key in signals.keys()]))
#cvarname = varnames[0] if 'error' not in varnames[0] else varnames[1]
#evarname = varnames[1] if 'error' not in varnames[0] else varnames[0]

for topo in topologies:
    if 'Int' in topo:
        topo_dict = {
            topo : {
                'pos' : {
                    'Central' : dense_lookup(*signals[('%s/%s_pos' % (topo, cvarname), 'dense_lookup')]),
                    'Error' : dense_lookup(*signals[('%s/%s_pos_error' % (topo, cvarname), 'dense_lookup')]),
                },
                'neg' : {
                    'Central' : dense_lookup(*signals[('%s/%s_neg' % (topo, cvarname), 'dense_lookup')]),
                    'Error' : dense_lookup(*signals[('%s/%s_neg_error' % (topo, cvarname), 'dense_lookup')]),
                },
            }
        }
    else:
        topo_dict = {
            topo : {
                'pos' : {
                    'Central' : dense_lookup(*signals[('%s/%s_pos' % (topo, cvarname), 'dense_lookup')]),
                    'Error' : dense_lookup(*signals[('%s/%s_pos_error' % (topo, cvarname), 'dense_lookup')]),
                },
            }
        }
    #set_trace()
    output_dict.update(topo_dict)

#set_trace()
outname = '%s.coffea' % args.outname if args.outname else 'signal_%s.coffea' % rfile.split('.')[0]
save(output_dict, outname)
print('%s written' % outname)
