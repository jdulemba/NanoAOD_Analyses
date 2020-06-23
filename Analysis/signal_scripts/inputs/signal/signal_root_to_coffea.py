from coffea.lookup_tools.root_converters import convert_histo_root_file
from coffea.lookup_tools.dense_lookup import dense_lookup
from pdb import set_trace
from coffea.util import save

rfile = 'signal_dists.root'

signals = convert_histo_root_file(rfile)

output_dict = {}

topologies = sorted(list(set([key[0].split('/')[0] for key in signals.keys()])))
varnames = list(set([key[0].split('/')[1] for key in signals.keys()]))
cvarname = varnames[0] if 'error' not in varnames[0] else varnames[1]
evarname = varnames[1] if 'error' not in varnames[0] else varnames[0]
for topo in topologies:
    topo_dict = {
        topo : {
            'Central' : dense_lookup(*signals[('%s/%s' % (topo, cvarname), 'dense_lookup')]),
            'Error' : dense_lookup(*signals[('%s/%s' % (topo, evarname), 'dense_lookup')]),
        }
    }
    central = dense_lookup(*signals[('%s/%s' % (topo, cvarname), 'dense_lookup')])
    error = dense_lookup(*signals[('%s/%s' % (topo, evarname), 'dense_lookup')])
    output_dict.update(topo_dict)

outname = 'signal_dists.coffea'
save(output_dict, outname)
print('%s written' % outname)
