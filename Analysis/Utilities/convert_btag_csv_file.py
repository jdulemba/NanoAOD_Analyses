import numpy
import sys
import warnings
from pdb import set_trace

# pt except for reshaping, then discriminant
btag_feval_dims = {"L": [1], "M": [1], "T": [1], "R": [2]}
#btag_feval_dims = {0: [1], 1: [1], 2: [1], 3: [2]}


def convert_btag_csv_file(csvFilePath):
    warnings.warn(
        "Auto-conversion of btag CSV files is deprecated. Try coffea.btag_tools.BTagScaleFactor!",
        FutureWarning,
    )

    fopen = open
    fmode = "rt"
    if ".gz" in csvFilePath:
        import gzip

        fopen = gzip.open
        fmode = (
            "r"
            if sys.platform.startswith("win") and sys.version_info.major < 3
            else fmode
        )
    btag_f = fopen(csvFilePath, fmode)
    nameandcols = btag_f.readline().split(";")
    btag_f.close()
    name = "btagsf"
    columns = None
    if len(nameandcols) == 2:
        name = nameandcols[0].strip()
        columns = nameandcols[1].strip()
    else:
        warnings.warn(
            "btagging SF file does not contain a name, using default!", RuntimeWarning
        )
        columns = nameandcols[0].strip()
    columns = [column.strip() for column in columns.split(",")]

    #set_trace()
    try:
        corrections = numpy.genfromtxt(
            csvFilePath,
            dtype=None,
            names=tuple(columns),
            converters={
                1: lambda s: s.strip(),
                2: lambda s: s.strip(),
                10: lambda s: s.strip(' "'),
            },
            delimiter=",",
            skip_header=1,
            encoding="ascii",
        )
    except:
        # needed to handle when 'pow' is used in formula
        import pandas
        df = pandas.read_csv(csvFilePath, skipinitialspace=True)
        df_numpy = df.to_numpy()
        corr_list = [tuple(df_numpy[idx]) for idx in range(len(df_numpy))]

            # hardcoded
        dtype_list = [('OperatingPoint', '<U1'), ('measurementType', '<U4'), ('sysType', '<U30'), ('jetFlavor', '<i8'), ('etaMin', '<f8'), ('etaMax', '<f8'), ('ptMin', '<f8'), ('ptMax', '<f8'), ('discrMin', '<f8'), ('discrMax', '<f8'), ('formula', '<U1000')]
        corrections = numpy.array(corr_list, dtype=dtype_list)
        
    #set_trace()
    all_names = corrections[[columns[i] for i in range(4)]]
    labels = numpy.unique(corrections[[columns[i] for i in range(4)]])
    wrapped_up = {}
    for label in labels:
        etaMins = numpy.unique(corrections[numpy.where(all_names == label)][columns[4]])
        etaMaxs = numpy.unique(corrections[numpy.where(all_names == label)][columns[5]])
        etaBins = numpy.union1d(etaMins, etaMaxs).astype(numpy.double)
        ptMins = numpy.unique(corrections[numpy.where(all_names == label)][columns[6]])
        ptMaxs = numpy.unique(corrections[numpy.where(all_names == label)][columns[7]])
        ptBins = numpy.union1d(ptMins, ptMaxs).astype(numpy.double)
        discrMins = numpy.unique(
            corrections[numpy.where(all_names == label)][columns[8]]
        )
        discrMaxs = numpy.unique(
            corrections[numpy.where(all_names == label)][columns[9]]
        )
        discrBins = numpy.union1d(discrMins, discrMaxs).astype(numpy.double)
        vals = numpy.zeros(
            shape=(len(discrBins) - 1, len(ptBins) - 1, len(etaBins) - 1),
            dtype=corrections.dtype[10],
        )
        for i, eta_bin in enumerate(etaBins[:-1]):
            for j, pt_bin in enumerate(ptBins[:-1]):
                for k, discr_bin in enumerate(discrBins[:-1]):
                    this_bin = numpy.where(
                        (all_names == label)
                        & (corrections[columns[4]] == eta_bin)
                        & (corrections[columns[6]] == pt_bin)
                        & (corrections[columns[8]] == discr_bin)
                    )
                    vals[k, j, i] = corrections[this_bin][columns[10]][0]
        label_decode = []
        for i in range(len(label)):
            label_decode.append(label[i])
            if isinstance(label_decode[i], bytes):
                label_decode[i] = label_decode[i].decode()
            else:
                label_decode[i] = str(label_decode[i])
        str_label = "_".join([name] + label_decode)
        feval_dim = btag_feval_dims[label[0]]
        wrapped_up[(str_label, "dense_evaluated_lookup")] = (
            vals,
            (etaBins, ptBins, discrBins),
            tuple(feval_dim),
        )
    return wrapped_up
