from pdb import set_trace
import fnmatch
import os
import matplotlib.pyplot as plt

base_jobid = os.environ['jobid']

#dataset_groups = {'2016APV' : {}, '2016' : {}, '2017' : {}, '2018' : {}}
dataset_groups = {
    year : {
        'EWK' : ['[WZ][WZ]', '[WZ]Jets', 'tt[WZ]*'],
        'ttJets_right' : ['ttJets*_right'],
        'ttJets_matchable' : ['ttJets*_matchable'],
        'ttJets_unmatchable' : ['ttJets*_unmatchable'],
        'ttJets_sl_tau' : ['ttJets*_sl_tau'],
        'ttJets_other' : ['ttJets*_other'],
        'ttJets' : ['ttJets', 'ttJets_PS'] if (base_jobid == 'NanoAODv6' and year == '2016') else ['ttJetsSL', 'ttJetsDiLep', 'ttJetsHad'],
        'ttJets_Sys' : ['ttJets_*DOWN', 'ttJets_*UP'] if (base_jobid == 'NanoAODv6' and year == '2016') else ['ttJetsSL_*DOWN', 'ttJetsDiLep_*DOWN', 'ttJetsHad_*DOWN', 'ttJetsSL_*UP', 'ttJetsDiLep_*UP', 'ttJetsHad_*UP'],
        'singlet' : ['single*'],
        'QCD' : ['QCD*'],
        'data' : ['data_Single*'],
    }
    for year in ['2016APV', '2016', '2017', '2018']
}
#dataset_groups['2016APV'] = dataset_groups['2016']
#dataset_groups['2017'] = {
#    'EWK' : ['[WZ][WZ]', '[WZ]Jets', 'tt[WZ]*'],
#    'ttJets_right' : ['ttJets*_right'],
#    'ttJets_matchable' : ['ttJets*_matchable'],
#    'ttJets_unmatchable' : ['ttJets*_unmatchable'],
#    'ttJets_sl_tau' : ['ttJets*_sl_tau'],
#    'ttJets_other' : ['ttJets*_other'],
#    'ttJets' : ['ttJetsSL', 'ttJetsDiLep', 'ttJetsHad'],#, 'ttJets_right', 'ttJets_matchable', 'ttJets_unmatchable', 'ttJets_other'],
#    'ttJets_Sys' : ['ttJetsSL_*DOWN', 'ttJetsDiLep_*DOWN', 'ttJetsHad_*DOWN', 'ttJetsSL_*UP', 'ttJetsDiLep_*UP', 'ttJetsHad_*UP'],
#    'singlet' : ['single*'],
#    'QCD' : ['QCD*'],
#    'data' : ['data_Single*'],
#}
#dataset_groups['2018'] = {
#    'EWK' : ['[WZ][WZ]', '[WZ]Jets', 'tt[WZ]*'],
#    'ttJets_right' : ['ttJets*_right'],
#    'ttJets_matchable' : ['ttJets*_matchable'],
#    'ttJets_unmatchable' : ['ttJets*_unmatchable'],
#    'ttJets_sl_tau' : ['ttJets*_sl_tau'],
#    'ttJets_other' : ['ttJets*_other'],
#    'ttJets' : ['ttJetsSL', 'ttJetsDiLep', 'ttJetsHad'],#, 'ttJets_right', 'ttJets_matchable', 'ttJets_unmatchable', 'ttJets_other'],
#    'ttJets_Sys' : ['ttJetsSL_*DOWN', 'ttJetsDiLep_*DOWN', 'ttJetsHad_*DOWN', 'ttJetsSL_*UP', 'ttJetsDiLep_*UP', 'ttJetsHad_*UP'],
#    'singlet' : ['single*'],
#    'QCD' : ['QCD*'],
#    'data' : ['data_Single*'],
#}

    ## dataset groupings for making templates to be used in fit
#template_groups = {'2016APV' : {}, '2016' : {}, '2017' : {}, '2018' : {}}
template_groups = {
    year: {
        'QCD' : ['QCD*'],
        'TT' : ['ttJets*_right', 'ttJets*_matchable', 'ttJets*_unmatchable', 'ttJets*_other'],
        'VV' : ['[WZ][WZ]'],
        'TTV' : ['tt[WZ]*'],
        'WJets' : ['WJets'],
        'ZJets' : ['ZJets'],
        'sChannel' : ['single*_schannel*'],
        'tChannel' : ['single*_tchannel*'],
        'tWChannel' : ['single*_tW*'],
        'data_obs' : ['data_Single*'],
    }
    for year in ['2016APV', '2016', '2017', '2018']
}
#template_groups['2017'] = {
#    'QCD' : ['QCD*'],
#    'TT' : ['ttJets*_right', 'ttJets*_matchable', 'ttJets*_unmatchable', 'ttJets*_other'],
#    'VV' : ['[WZ][WZ]'],
#    'TTV' : ['tt[WZ]*'],
#    'WJets' : ['WJets'],
#    'ZJets' : ['ZJets'],
#    'sChannel' : ['single*_schannel*'],
#    'tChannel' : ['single*_tchannel*'],
#    'tWChannel' : ['single*_tW*'],
#    'data_obs' : ['data_Single*'],
#}
#template_groups['2018'] = {
#    'QCD' : ['QCD*'],
#    'TT' : ['ttJets*_right', 'ttJets*_matchable', 'ttJets*_unmatchable', 'ttJets*_other'],
#    'VV' : ['[WZ][WZ]'],
#    'TTV' : ['tt[WZ]*'],
#    'WJets' : ['WJets'],
#    'ZJets' : ['ZJets'],
#    'sChannel' : ['single*_schannel*'],
#    'tChannel' : ['single*_tchannel*'],
#    'tWChannel' : ['single*_tW*'],
#    'data_obs' : ['data_Single*'],
#}

# create groups for signal samples and add them to dataset/template_groups
signal_groups = {
    'AtoTTJetsDiLep_M400_W10_Int' : ['AtoTTJetsDiLep_M400_W10_Int_*'],
    'HtoTTJetsDiLep_M750_W2p5_Int': ['HtoTTJetsDiLep_M750_W2p5_Int_*'],
    'AtoTTJetsDiLep_M400_W10_Res' : ['AtoTTJetsDiLep_M400_W10_Res'],
    'HtoTTJetsDiLep_M750_W2p5_Res': ['HtoTTJetsDiLep_M750_W2p5_Res'],
}
for group in dataset_groups.values():
    group.update(signal_groups)
for group in template_groups.values():
    group.update(signal_groups)


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

def get_label(sample, styles):
    best_pattern = ''
    for pattern, style_dict in styles.items():
        if fnmatch.fnmatch(sample, pattern):
            if len(pattern) > len(best_pattern):
                best_pattern = pattern
    if best_pattern:
        return styles[best_pattern]['name']
    else:
        return sample

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

def make_dataset_groups(lepton, year, samples=[], gdict='dataset'):
    proj_dir = os.environ['PROJECT_DIR']
    jobid = os.environ['jobid']

    if not samples:
            ## get file that has names for all datasets to use
        fpath = '/'.join([proj_dir, 'inputs', '%s_%s' % (year, jobid), 'analyzer_inputs.txt'])
        if not os.path.isfile(fpath):
            raise IOError("File %s not found" % fpath)

        txt_file = open(fpath, 'r')
        samples = [sample.strip('\n') for sample in txt_file if not sample.startswith('#')]
        samples = [sample for sample in samples if not (sample.startswith('data') and lepton not in sample)] # get rid of data samples that don't correspond to lepton chosen

    groupings = {}
    if gdict == 'dataset':
        groups_dict = dataset_groups 
    elif gdict == 'templates':
        groups_dict = template_groups
    else:
        raise ValueError("Can only choose between 'dataset' and 'templates' as inputs for gdict")

    for group_name, patterns in groups_dict[year].items():
        flist = []
        for sample in samples:
            if any([fnmatch.fnmatch(sample, pattern) for pattern in patterns]):
                flist.append(sample)
        if flist:
            groupings[group_name] = flist

    return groupings


def add_coffea_files(input_files):
    from coffea.util import load

    output_acc = load(input_files[0])
    for idx in range(1, len(input_files)):
        try:
            output_acc.add(load(input_files[idx]))
        except:
            raise ValueError("File %s (number %i) could not be added" % (input_files[idx], idx))

    return output_acc

def save_accumulator(accumulator, output_fname):
    from coffea.util import save
    save(accumulator, output_fname)
    print('%s written' % output_fname)


def print_table(lines, filename='', separate_header=True, header_line=0, print_output=False):
    ''' Prints a formatted table given a 2 dimensional array '''
    #Count the column width
    widths = []
    for line in lines:
        for i, size in enumerate( [ len(x) for x in line ] ):
            while i >= len(widths):
                widths.append(0)
            if size > widths[i]:
                widths[i] = size
    
    #Generate the format string to pad the columns
    print_string = ""
    for i, width in enumerate(widths):
        print_string += "{" + str(i) + ":" + str(width) + "} | "
    if (len(print_string) == 0):
        return
    print_string = print_string[:-3]
    
    if print_output: ## print data
        for i,line in enumerate(lines):
            print(print_string.format(*line))
            if (i == header_line and separate_header):
                print("-"*(sum(widths)+3*(len(widths)-1)))

    if filename != '': # save to filename    
        with open(filename, 'w') as f:
            #Write the data into filename
            for i,line in enumerate(lines):
                print(print_string.format(*line), file=f)
                if (i == header_line and separate_header):
                    print("-"*(sum(widths)+3*(len(widths)-1)), file=f)

