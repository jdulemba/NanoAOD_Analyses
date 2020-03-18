import os
from pdb import set_trace
import fnmatch
import Utilities.prettyjson as prettyjson

def get_sample_list(indir, sample=None, text_file=None):

    if sample and text_file:
        raise IOError("Only sample OR text file with samples to use should be input")
    if not (sample or text_file):
        raise IOError("Sample OR text file with samples to use must be input")

        ## get samples to use
    if sample:
            ## sample specified
        if not os.path.isfile('%s/%s.txt' % (indir, sample)):
            raise IOError("File with samples %s.txt not found" % sample)
    
        samples = ['%s/%s.txt' % (indir, sample)]
    
    else:
            ## get file that has names for all datasets to use
        fpath = '/'.join([indir, text_file])
        if not os.path.isfile(fpath):
            raise IOError("File %s not found" % text_file)
    
        txt_file = open(fpath, 'r')
        samples = ['%s/%s.txt' % (indir, sample.strip('\n')) for sample in txt_file if not sample.startswith('#')]
        if not samples:
            raise IOError("No samples found as inputs")

    return samples

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def get_file_range(lst, n):
    '''
    Input: list of integers to be broken up into groups of size n
    '''
    file_chunks = list(chunks(lst, n))
    franges = ['%i:%i' % (chunk[0], chunk[-1]) for chunk in file_chunks]
    return franges

def get_file_splitting(sample):
    splittings = prettyjson.loads(open('%s/Run_Jobs/splittings.json' % os.environ['PROJECT_DIR']).read())
    if '*' in splittings.keys():
        f_split = splittings['*']
    elif b'*' in splittings.keys():
        f_split = splittings[b'*']
    for pattern, splitting in splittings.items():
        if pattern == b'*': continue
        if fnmatch.fnmatch(sample, pattern):
            f_split = splitting
    return f_split
