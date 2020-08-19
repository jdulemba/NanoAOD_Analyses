from coffea.util import load
from pdb import set_trace
import os
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
from coffea import hist
import numpy as np
import fnmatch
import Utilities.Plotter as Plotter
import uproot
    
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016', '2017', '2018'], help='What year is the ntuple from.')
parser.add_argument('--maskData', action='store_false', help='Mask templates for data, default is True.')

args = parser.parse_args()

def get_templates(tmp_rname):
    '''
    Function that writes linearized mtt vs costheta distributions to root file.
    '''
    ## variables that only need to be defined/evaluated once
    hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])

        ## get data lumi and scale MC by lumi
    data_lumi_year = prettyjson.loads(open('%s/inputs/lumis_data.json' % proj_dir).read())[args.year]

        # get correct hist and rebin
    hname_to_use = 'mtt_vs_tlep_ctstar_abs'
    if hname_to_use not in hdict.keys():
        raise ValueError("%s not found in file" % hname_to_use)
    xrebinning, yrebinning = linearize_binning
    histo = hdict[hname_to_use] # process, sys, jmult, leptype, btag, lepcat
    
    xaxis_name = histo.dense_axes()[0].name
    yaxis_name = histo.dense_axes()[1].name
        ## rebin x axis
    if isinstance(xrebinning, np.ndarray):
        new_xbins = hist.Bin(xaxis_name, xaxis_name, xrebinning)
    elif isinstance(xrebinning, float) or isinstance(xrebinning, int):
        new_xbins = xrebinning
    histo = histo.rebin(xaxis_name, new_xbins)
        ## rebin y axis
    if isinstance(yrebinning, np.ndarray):
        new_ybins = hist.Bin(yaxis_name, yaxis_name, yrebinning)
    elif isinstance(yrebinning, float) or isinstance(yrebinning, int):
        new_ybins = yrebinning
    rebin_histo = histo.rebin(yaxis_name, new_ybins)
    
    nbins = (len(xrebinning)-1)*(len(yrebinning)-1)
    
        
        ## scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
    ttJets_permcats = ['*right', '*matchable', '*unmatchable', '*other']
    names = [dataset for dataset in sorted(set([key[0] for key in hdict[hname_to_use].values().keys()]))] # get dataset names in hists
    ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...
    
        ## make groups based on process
    process = hist.Cat("process", "Process", sorting='placement')
    process_cat = "dataset"

        # need to save coffea hist objects to file so they can be opened by uproot in the proper format
    upfout = uproot.recreate(tmp_rname, compression=uproot.ZLIB(4)) if os.path.isfile(tmp_rname) else uproot.create(tmp_rname)

    for lep in ['Muon', 'Electron']:
        lepdir = 'mujets' if lep == 'Muon' else 'ejets'
    
        ## make groups based on process
        process_groups = plt_tools.make_dataset_groups(lep, args.year, samples=names, gdict='templates')
        
        lumi_correction = load('%s/Corrections/%s/MC_LumiWeights_IgnoreSigEvts.coffea' % (proj_dir, jobid))[args.year]['%ss' % lep]
                # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
        if len(ttJets_cats) > 0:
            for tt_cat in ttJets_cats:
                ttJets_lumi_topo = '_'.join(tt_cat.split('_')[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
                ttJets_eff_lumi = lumi_correction[ttJets_lumi_topo]
                lumi_correction.update({tt_cat: ttJets_eff_lumi})
    
        histo = rebin_histo.copy()
        histo.scale(lumi_correction, axis='dataset')
        histo = histo.group(process_cat, process, process_groups)[:, :, :, lep, :, :].integrate('leptype')

        #for jmult in ['3Jets']:        
        for jmult in sorted(set([key[2] for key in histo.values().keys()])):
            iso_sb    = Plotter.linearize_hist(histo[:, 'nosys', jmult, 'btagPass', 'Loose'].integrate('sys').integrate('jmult').integrate('lepcat').integrate('btag'))
            btag_sb   = Plotter.linearize_hist(histo[:, 'nosys', jmult, 'btagFail', 'Tight'].integrate('sys').integrate('jmult').integrate('lepcat').integrate('btag'))
            double_sb = Plotter.linearize_hist(histo[:, 'nosys', jmult, 'btagFail', 'Loose'].integrate('sys').integrate('jmult').integrate('lepcat').integrate('btag'))
            sig_histo = Plotter.linearize_hist(histo[:, :, jmult, 'btagPass', 'Tight'].integrate('jmult').integrate('lepcat').integrate('btag'))
        
            for sys in systematics.keys():
                if sys not in histo.axis('sys')._sorted:
                    print('\n\n   Systematic %s not available, skipping\n\n' % sys)
                    continue
                sysname, onlyTT = systematics[sys]
                if 'Lep_RECO' in sys: sysname = sysname.replace('LEP', lepdir[0])
        
                qcd_est_histo = Plotter.QCD_Est(sig_reg=sig_histo, iso_sb=iso_sb, btag_sb=btag_sb, double_sb=double_sb, norm_type='Sideband', shape_region='BTAG', norm_region='BTAG', sys=sys)
        
                for proc in sorted(set([key[0] for key in qcd_est_histo.values().keys()])):
                    if (proc != 'TT') and onlyTT: continue
                    if (proc == 'data_obs') and not (sys == 'nosys'): continue
                    name = proc+lepdir if proc == 'QCD' else proc
                    print(lep, jmult, sys, name)
                    #set_trace()
                    outhname = '_'.join([jmult, lepdir, name]) if sys == 'nosys' else '_'.join([jmult, lepdir, name, sysname])
                    template_histo = qcd_est_histo[proc].integrate('process')
                    #if proc == 'data_obs': set_trace()
                    upfout[outhname] = hist.export1d(template_histo)
    
    upfout.close()
    print('%s written' % tmp_rname)


def write_correct_template_format(in_fname):
    '''
    Opens temporary root file where template distributions are and then saves them with the correct structure/naming
    '''

    from rootpy.io import root_open
    
    rfile = root_open(in_fname) if in_fname.endswith('.root') else root_open('%s.root' % in_fname)
    
    mu_3j_keys = [key.name for key in rfile.keys() if '3Jets_mujets' in key.name]
    el_3j_keys = [key.name for key in rfile.keys() if '3Jets_ejets' in key.name]
    mu_4pj_keys = [key.name for key in rfile.keys() if '4PJets_mujets' in key.name]
    el_4pj_keys = [key.name for key in rfile.keys() if '4PJets_ejets' in key.name]
    
    #set_trace()
    fname_3j = '%s/templates_lj_3Jets_bkg_%s_QCD_Est_%s.root' % (outdir, jobid, args.year)
    with root_open(fname_3j, 'w') as rout:
        mu_dir = rout.mkdir('mujets')
        mu_dir.cd()

        for key in mu_3j_keys:
            hname = key.split('3Jets_mujets_')[-1]
            histo = rfile.Get(key)
            if (hname == 'data_obs') and (args.maskData):
                histo.Reset()
            mu_dir.WriteTObject(histo, hname)
    
        el_dir = rout.mkdir('ejets')
        el_dir.cd()
    
        for key in el_3j_keys:
            hname = key.split('3Jets_ejets_')[-1]
            histo = rfile.Get(key)
            if (hname == 'data_obs') and (args.maskData):
                histo.Reset()
            el_dir.WriteTObject(histo, hname)

    print('%s written' % fname_3j)

    fname_4pj = '%s/templates_lj_4PJets_bkg_%s_QCD_Est_%s.root' % (outdir, jobid, args.year)
    with root_open(fname_4pj, 'w') as rout:
        mu_dir = rout.mkdir('mujets')
        mu_dir.cd()
        for key in mu_4pj_keys:
            hname = key.split('4PJets_mujets_')[-1]
            histo = rfile.Get(key)
            if (hname == 'data_obs') and (args.maskData):
                histo.Reset()
            mu_dir.WriteTObject(histo, hname)
    
        el_dir = rout.mkdir('ejets')
        el_dir.cd()
        for key in el_4pj_keys:
            hname = key.split('4PJets_ejets_')[-1]
            histo = rfile.Get(key)
            if (hname == 'data_obs') and (args.maskData):
                histo.Reset()
            el_dir.WriteTObject(histo, hname)

    print('%s written' % fname_4pj)



if __name__ == '__main__':
    proj_dir = os.environ['PROJECT_DIR']
    jobid = os.environ['jobid']
    analyzer = 'htt_btag_iso_cut'
    
    input_dir = '/'.join([proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer])
    f_ext = 'TOT.coffea'
    
    fnames = sorted(['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])
    
    outdir = '/'.join([proj_dir, 'plots', '%s_%s' % (args.year, jobid), analyzer, 'Templates'])
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    
        # (name, only TT)
    systematics = {
        'nosys' : ('', False),
        'JES_UP' : ('CMS_scale_j_13TeVUp', False),
        'JES_DW' : ('CMS_scale_j_13TeVDown', False),
        'JER_UP' : ('CMS_res_j_13TeVUp', False),
        'JER_DW' : ('CMS_res_j_13TeVDown', False),
        'btag_bc_UP' : ('CMS_eff_b_13TeVUp', False),
        'btag_bc_DW' : ('CMS_eff_b_13TeVDown', False),
        'btag_l_UP' : ('CMS_fake_b_13TeVUp', False),
        'btag_l_DW' : ('CMS_fake_b_13TeVDown', False),
        'Lep_RECOUp' : ('CMS_eff_LEPUp', False),
        'Lep_RECODown' : ('CMS_eff_LEPDown', False),
        'PileupUp' : ('CMS_pileupUp', False),
        'PileupDown' : ('CMS_pileupDown', False),
        'MET_UP' : ('CMS_METunclustered_13TeVUp', False),
        'MET_DW' : ('CMS_METunclustered_13TeVDown', False),
        'RENORMUp' : ('QCDscaleMERenorm_TTUp', True),
        'RENORMDown' : ('QCDscaleMERenorm_TTDown', True),
        'FACTORUp' : ('QCDscaleMEFactor_TTUp', True),
        'FACTORDown' : ('QCDscaleMEFactor_TTDown', True),
        'RENORM_FACTOR_SAMEUp' : ('QCDscaleMERenormFactor_TTUp', True),
        'RENORM_FACTOR_SAMEDown' : ('QCDscaleMERenormFactor_TTDown', True),
        'hdampUP' : ('Hdamp_TTUp', True),
        'hdampDOWN' : ('Hdamp_TTDown', True),
        'mtopUP' : ('TMassUp', True),
        'mtopDOWN' : ('TMassDown', True),
    }
    
    
    linearize_binning = (
        np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 625.0, 650.0, 675.0, 700.0, 730.0, 760.0, 800.0, 850.0, 900.0, 1000.0, 1200.0]),
        #np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 650., 700.0, 750.0, 800.0, 850.0, 900.0, 2000.0]),
        np.array([0.0, 0.4, 0.6, 0.75, 0.9, 1.0])
    )

    temp_rname = 'tmp.root'
    try:
        get_templates(temp_rname)
        #set_trace()
        write_correct_template_format(temp_rname)
        os.system('rm %s' % temp_rname)
        print('%s deleted' % temp_rname)
        
    except:
        print('Could not write templates to file')
