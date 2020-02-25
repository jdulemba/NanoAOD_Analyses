from pdb import set_trace
import fnmatch


dataset_groups = {
    'EWK' : ['[WZ][WZ]', '[WZ]Jets', 'tt[WZ]*'],
    'ttJets' : ['ttJets', 'ttJets_PS'],#, 'ttJets_right', 'ttJets_matchable', 'ttJets_unmatchable', 'ttJets_other'],
    'ttJets_Sys' : ['ttJets_*DOWN', 'ttJets_*UP'],
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
        return None

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

