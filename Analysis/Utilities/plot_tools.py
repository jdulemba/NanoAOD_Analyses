from pdb import set_trace
import fnmatch
import os

dataset_groups = { '2016' : {}, '2017' : {}, '2018' : {}}
dataset_groups['2016'] = {
    'EWK' : ['[WZ][WZ]', '[WZ]Jets', 'tt[WZ]*'],
    'ttJets' : ['ttJets', 'ttJets_PS'],#, 'ttJets_right', 'ttJets_matchable', 'ttJets_unmatchable', 'ttJets_other'],
    'ttJets_Sys' : ['ttJets_*DOWN', 'ttJets_*UP'],
    'singlet' : ['single*'],
    'QCD' : ['QCD*'],
    'data' : ['data_Single*'],
}
dataset_groups['2017'] = {
    'EWK' : ['[WZ][WZ]', '[WZ]Jets', 'tt[WZ]*'],
    'ttJets' : ['ttJetsSL', 'ttJetsDiLep', 'ttJetsHad'],#, 'ttJets_right', 'ttJets_matchable', 'ttJets_unmatchable', 'ttJets_other'],
    'ttJets_Sys' : ['ttJetsSL_*', 'ttJetsDiLep_*', 'ttJetsHad_*'],
    'singlet' : ['single*'],
    'QCD' : ['QCD*'],
    'data' : ['data_Single*'],
}
dataset_groups['2018'] = {
    'EWK' : ['[WZ][WZ]', '[WZ]Jets', 'tt[WZ]*'],
    'ttJets' : ['ttJetsSL', 'ttJetsDiLep', 'ttJetsHad'],#, 'ttJets_right', 'ttJets_matchable', 'ttJets_unmatchable', 'ttJets_other'],
    'ttJets_Sys' : ['ttJetsSL_*', 'ttJetsDiLep_*', 'ttJetsHad_*'],
    'singlet' : ['single*'],
    'QCD' : ['QCD*'],
    'data' : ['data_Single*'],
}

def get_styles(sample, styles):
    best_pattern = ''
    for pattern, style_dict in styles.items():
        if fnmatch.fnmatch(sample, pattern):
            if len(pattern) > len(best_pattern):
                best_pattern = pattern
    if best_pattern:
        return styles[best_pattern]['facecolor'], styles[best_pattern]['name']
    else:
        return 'r', sample

def get_group(sample, styles=dataset_groups):
    best_pattern = ''
    for group_name, patterns in styles.items():
        if any([fnmatch.fnmatch(sample, pattern) for pattern in patterns]):
            best_pattern = group_name

    if best_pattern:
        return best_pattern
    else:
        print("Pattern not found for %s" % sample)
        return sample

def make_dataset_groups(lepton, year):
    proj_dir = os.environ['PROJECT_DIR']
    jobid = os.environ['jobid']

        ## get file that has names for all datasets to use
    fpath = '/'.join([proj_dir, 'inputs', '%s_%s' % (year, jobid), 'analyzer_inputs.txt'])
    if not os.path.isfile(fpath):
        raise IOError("File %s not found" % fpath)

    txt_file = open(fpath, 'r')
    samples = [sample.strip('\n') for sample in txt_file if not sample.startswith('#')]
    samples = [sample for sample in samples if not (sample.startswith('data') and lepton not in sample)] # get rid of data samples that don't correspond to lepton chosen

    groupings = {}
    for group_name, patterns in dataset_groups[year].items():
        flist = []
        for sample in samples:
            if any([fnmatch.fnmatch(sample, pattern) for pattern in patterns]):
                flist.append(sample)
        if flist:
            groupings[group_name] = flist

    return groupings

#groupings = make_dataset_groups('Muon', '2016')


def add_coffea_files(input_files):
    from coffea.util import load
    import collections
    input_accs = [load(fname) for fname in input_files]

    output_acc = collections.Counter()
    for acc in input_accs:
        output_acc.update(acc)

    return output_acc

def save_accumulator(accumulator, output_fname):
    from coffea.util import save
    save(accumulator, output_fname)
    print('%s written' % output_fname)

