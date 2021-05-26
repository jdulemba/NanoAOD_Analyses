import os
from pdb import set_trace
import Utilities.prettyjson as prettyjson

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('--type', default='data', nargs='?', choices=['SM', 'signal', 'data', 'all'], help='specify which sample type to print')
args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']

    # only print signal MC info
if (args.type == 'signal') or (args.type == 'all'):
    import itertools
    from string import Template
    name2val = lambda x: float(x.replace('pc','').replace('p', '.'))

    samples_list = prettyjson.loads(open('%s/inputs/samples_2016.json' % proj_dir).read())
    
    samples_dict = {}
    for sample in samples_list:
        name = sample.pop('name')
        samples_dict[name] = sample

        ## make table for A
    A_output = "\multirow{2}{*}{Parity} & \multirow{2}{*}{$\mathsf{m_{A}}$ [GeV]} & \multirow{2}{*}{$\Gamma_{\mathsf{A}}$ [\% $\mathsf{m_{A}}$]} & \multicolumn{2}{c |}{LO $\sigma$ [pb]} & \multirow{2}{*}{$\mathsf{k_{R}}$} \\\ \n"
    A_output += " & & & Resonance & Interference & \\\ \n\hline \n"
    H_output = "\multirow{2}{*}{Parity} & \multirow{2}{*}{$\mathsf{m_{H}}$ [GeV]} & \multirow{2}{*}{$\Gamma_{\mathsf{H}}$ [\% $\mathsf{m_{H}}$]} & \multicolumn{2}{c |}{LO $\sigma$ [pb]} & \multirow{2}{*}{$\mathsf{k_{R}}$} \\\ \n"
    H_output += " & & & Resonance & Interference & \\\ \n\hline \n"
    for sig_point in itertools.product(['M400','M500', 'M600', 'M750'], ['W2p5', 'W5', 'W10', 'W25']):
        signal = '_'.join(sig_point)

        mass = sig_point[0][1:]
        width = name2val(sig_point[1][1:])

        A_res_xsec = samples_dict['AtoTT_%s_Res' % signal]['xsection']
        A_int_xsec = samples_dict['AtoTT_%s_Int' % signal]['xsection']
        H_res_xsec = samples_dict['HtoTT_%s_Res' % signal]['xsection']
        H_int_xsec = samples_dict['HtoTT_%s_Int' % signal]['xsection']

        A_template = Template("\multirow{16}{*}{A} & \multirow{4}{*}{$MASS} & $WIDTH & $XSEC_A_RES & $XSEC_A_INT & $KFACTOR \\\\ \n")
        A_output += A_template.substitute(MASS=mass, WIDTH=width, XSEC_A_RES="{:.3f}".format(A_res_xsec), XSEC_A_INT="{:.3f}".format(A_int_xsec), KFACTOR='---')
        H_template = Template("\multirow{16}{*}{H} & \multirow{4}{*}{$MASS} & $WIDTH & $XSEC_H_RES & $XSEC_H_INT & $KFACTOR \\\\ \n")
        H_output += H_template.substitute(MASS=mass, WIDTH=width, XSEC_H_RES="{:.3f}".format(H_res_xsec), XSEC_H_INT="{:.3f}".format(H_int_xsec), KFACTOR='---')

    print(A_output)
    print(H_output)
    sig_datasets_info = open('%s/AN_Scripts/sig_datasets.txt' % proj_dir, 'w')
    sig_datasets_info.write(A_output)
    sig_datasets_info.write(H_output)
    sig_datasets_info.close()
    print('%s/AN_Scripts/sig_datasets.txt written' % proj_dir)


lumi_list = prettyjson.loads(open('%s/inputs/lumis_data.json' % proj_dir).read())

years = ['2016', '2017', '2018']
for year in years:
    samples_list = prettyjson.loads(open('%s/inputs/samples_%s.json' % (proj_dir, year)).read())
    
    samples_dict = {}
    for sample in samples_list:
        name = sample.pop('name')
        samples_dict[name] = sample

    if (args.type == 'data') or (args.type == 'all'):
        lumiMasks = {}    
        data_output = "Data set & Integrated luminosity ($\\fbinv$)\\\\ \n\hline \n"
        data_output += "{YEAR} & {MU_LUMI:.2f} ({EL_LUMI:.2f})\\\\ \n".format(YEAR=year, MU_LUMI=lumi_list[year]['Muons']/1000, EL_LUMI=lumi_list[year]['Electrons']/1000)
    
    if (args.type == 'SM') or (args.type == 'all'):
        sm_output = "Data set & $\sigma$ [pb] & Weights\\\\ \n\hline \n"


    for sample in samples_dict.keys():
            # only print SM MC info
        if (args.type == 'SM') or (args.type == 'all'):
            if (sample.startswith('AtoTT')) or (sample.startswith('HtoTT')): continue
            if 'data' in sample: continue
            dbs_name = samples_dict[sample]['DBSName'] # get full DBS name
            name_to_use = dbs_name.split('/')[1] if dbs_name.startswith('/') else dbs_name.split('/')[0]
            xsec = samples_dict[sample]['xsection']
            nwts = '-' if not os.path.isfile(os.path.join(proj_dir, 'inputs', '%s_%s' % (year, jobid), '%s.meta.json' % sample)) else prettyjson.loads(open('%s/inputs/%s_%s/%s.meta.json' % (proj_dir, year, jobid, sample)).read())["nWeightedEvts"]
        
            sm_output += "{SAMPLE} & {XSEC} & {WTS} \\\\ \n".format(SAMPLE=name_to_use.replace('_', '\_'), XSEC=xsec, WTS=nwts)


            # only print data info
        if (args.type == 'data') or (args.type == 'all'):
            if not 'data' in sample: continue
            lepton = sample.split('data_Single')[-1].split('_')[0]+'s'
            period = sample.split(year)[-1]
            if period not in lumi_list['Ind'][year][lepton].keys():
                print('%s not found in list of lumis used. Make sure this is correct' % sample)
                continue
            lumi = lumi_list['Ind'][year][lepton][period]/1000
            dbs_name = samples_dict[sample]['DBSName'] # get full DBS name
            dbs_name = '/'.join(dbs_name.split('/')[:-1]) # remove 'NANOAOD' part of dbs name
            data_output += "{SAMPLE} & {LUMI:.2f} \\\\ \n".format(SAMPLE=dbs_name.replace('_', '\_'), LUMI=lumi)

                # add lumimask used for datset
            lumiMasks[sample] = samples_dict[sample]['lumimask']
    
    if (args.type == 'data') or (args.type == 'all'):
            # add lumimask to output
        if len(sorted(set(lumiMasks.values()))) == 1:
            data_output += "\n\n%s\n\n" % sorted(set(lumiMasks.values()))[0]

        print(data_output)
        data_datasets_lumi = open('%s/AN_Scripts/%s_data_datasets.txt' % (proj_dir, year), 'w')
        data_datasets_lumi.write(data_output)
        data_datasets_lumi.close()
        print('%s/AN_Scripts/%s_data_datasets.txt written' % (proj_dir, year))

    if (args.type == 'SM') or (args.type == 'all'):
        SM_datasets_info = open('%s/AN_Scripts/%s_SM_datasets.txt' % (proj_dir, year), 'w')
        SM_datasets_info.write(sm_output)
        SM_datasets_info.close()
        print('%s/AN_Scripts/%s_SM_datasets.txt written' % (proj_dir, year))
