from pdb import set_trace
import subprocess as sub

def query(query_str, verbose=False):
    'simple query function to interface with DAS'
    if verbose:
        print('querying DAS with: "%s"' % query_str)
    cmd=[
        'dasgoclient',
        '--query=%s' % query_str,
        '--limit=0',
        '--format=plain'
        ]
    p = sub.Popen(cmd, stdout=sub.PIPE,stderr=sub.PIPE)
    output, errors = p.communicate()
    if errors:
        raise RuntimeError('Das query crashed with error: %s' % errors)
    #-1 works both when getting dataset from files and files from datasets,
    #not checked on everything
    output = output.decode() # convert from bytes to string
    return output.split('\n')[:-1]
