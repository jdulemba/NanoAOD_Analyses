import os
from pdb import set_trace
import fnmatch
import Run_Jobs.splittings as splitting_dicts

def get_sample_list(indir, sample=None, text_file=None):

    if sample and text_file:
        raise IOError("Only sample OR text file with samples to use should be input")
    if not (sample or text_file):
        raise IOError("Sample OR text file with samples to use must be input")

        ## get samples to use
    if sample:
        ## match sample to any available in 'indir'
        fname_opts = [fname.split(".")[0] for fname in os.listdir(indir) if fname.endswith(".txt")]
        samples = ["%s/%s.txt" % (indir, name) for name in fname_opts if fnmatch.fnmatch(name, sample)]

    else:
            ## get file that has names for all datasets to use
        fpath = os.path.join(indir, text_file)
        if not os.path.isfile(fpath):
            raise IOError(f"File {fpath} not found")
    
        txt_file = open(fpath, "r")
        samples = ["%s/%s.txt" % (indir, sample.strip("\n")) for sample in txt_file if not sample.startswith("#")]

    if not samples:
        raise ValueError("No samples found to match %s" % sample)

    return samples

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def get_file_range(lst, n):
    """
    Input: list of integers to be broken up into groups of size n
    """
    file_chunks = list(chunks(lst, n))
    franges = ["%i:%i" % (chunk[0], chunk[-1]) for chunk in file_chunks]
    return franges

def get_file_splitting(sample, analyzer):
    if (analyzer == "meta_info") or ("check_" in analyzer):
        splittings = splitting_dicts.meta_dict
    elif analyzer == "htt_flav_effs":
        splittings = splitting_dicts.hflav_dict
    elif (analyzer == "ttbar_alpha_reco") or (analyzer == "ttbar_post_alpha_reco"):
        splittings = splitting_dicts.alpha_dict
    elif analyzer == "permProbComputer":
        splittings = splitting_dicts.perm_dict
    else:
        splittings = splitting_dicts.other_dict

    if "*" in splittings.keys():
        f_split = splittings["*"]
    elif b"*" in splittings.keys():
        f_split = splittings[b"*"]
    for pattern, splitting in splittings.items():
        if pattern == b"*": continue
        if fnmatch.fnmatch(sample, pattern):
            f_split = splitting
    return f_split
