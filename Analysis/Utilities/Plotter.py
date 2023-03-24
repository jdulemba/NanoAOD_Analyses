from coffea.hist import plot
# matplotlib
import matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
from Utilities import styles
import Utilities.plot_tools as plt_tools
import re
from pdb import set_trace
import numpy as np
import Utilities.systematics as systematics
from scipy import stats

## make data and mc categories for data/MC plotting
mc_samples = re.compile("(?!data*)")
data_samples = re.compile("(data*)")
qcd_samples = re.compile("(QCD*)")
prompt_mc_mask = re.compile(r"(?!(?:%s))" % "|".join(["data*", "QCD*"]))
nonTT_mc_mask = re.compile(r"(?!(?:%s))" % "|".join(["data*", "ttJets*", "TT$"])) # "TT$" is for templates
signal_samples = re.compile(r"((?:%s))" % "|".join(["AtoTT*", "HtoTT*"]))
nonsignal_samples = re.compile(r"(?!(?:%s))" % "|".join(["AtoTT*", "HtoTT*"]))
top_samples = re.compile(r"((?:%s))" % "|".join(["ttJets*", "singlet*"]))
bkg_samples = re.compile(r"((?:%s))" % "|".join(["QCD*", "EWK*"]))
#bkg_mask = re.compile("(BKG*)")
bkg_mask = re.compile("(EWQCD*)")
template_top_samples = re.compile(r"((?:%s))" % "|".join(["TT$", "TQ$", "TW$", "TB$"]))

hstyles = styles.styles
stack_fill_opts = {"alpha" : 0.8, "edgecolor" :( 0, 0, 0, .5), "linewidth" : 1.0}
stack_error_opts = {"edgecolor":(0,0,0,.5)}

def plot_stack1d(ax, rax, hdict, xlabel="", ylabel="", sys="nosys", xlimits=None, ylimits=None, **mc_opts):

    mcorder = mc_opts.get("mcorder")
    maskData = mc_opts.get("maskData", False)
    overflow = mc_opts.get("overflow", "none")
    logy = mc_opts.get("logy", False)
    leg_ncols = mc_opts.get("leg_ncols", 3)

    #set_trace()
    if hdict.sparse_dim() > 1:
        mc_dict = hdict[:, sys].integrate("sys")
        mc_dict = mc_dict[mc_samples]
        data_dict = hdict[:, "nosys"].integrate("sys")
        data_dict = data_dict[data_samples]
    
    else:
        mc_dict = hdict[mc_samples]
        data_dict = hdict[data_samples]

    #set_trace()    
        ## plot MC and data
    plot.plot1d(mc_dict,
        overlay=mc_dict.axes()[0].name,
        ax=ax,
        clear=False,
        stack=True,
        line_opts=None,
        fill_opts=stack_fill_opts,
        error_opts=stack_error_opts,
        order=mcorder,
        overflow=overflow,
        overlay_overflow=overflow,
    )
    if not maskData:
        plot.plot1d(data_dict,
            overlay=data_dict.axes()[0].name,
            ax=ax,
            clear=False,
            error_opts=hstyles["data_err_opts"],
            overflow=overflow,
            overlay_overflow=overflow,
        )
    if logy:
        ax.set_yscale("log")
    ax.autoscale()
    ax.set_ylim(0, ax.get_ylim()[1]*1.3)
    #ax.set_ylim(0, None)
    ax.set_xlabel(None)
    ax.set_xlim(xlimits)

        ## set legend and corresponding colors
    handles, labels = ax.get_legend_handles_labels()
    #set_trace()    
    for idx, sample in enumerate(labels):
        if sample == "data" or sample == "Observed": continue
        if isinstance(handles[idx], matplotlib.lines.Line2D): continue
        facecolor, legname = plt_tools.get_styles(sample, hstyles)
        handles[idx].set_facecolor(facecolor)
        labels[idx] = legname
    # call ax.legend() with the new values
    ax.legend(handles,labels, loc="upper right", title=mc_opts["legend_title"], ncol=leg_ncols) if "legend_title" in mc_opts.keys() else ax.legend(handles,labels, loc="upper right", ncol=leg_ncols)
    #ax.legend(handles,labels, loc="upper right", title=mc_opts["legend_title"]) if "legend_title" in mc_opts.keys() else ax.legend(handles,labels, loc="upper right")
    #set_trace()
    
        ## plot data/MC ratio
    if maskData:
        rax.axhspan(0.5, 1.5, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
    else:
        plot.plotratio(data_dict.sum(data_dict.axes()[0].name), mc_dict.sum(mc_dict.axes()[0].name),
            ax=rax,
            error_opts=hstyles["data_err_opts"],
            denom_fill_opts={},
            guide_opts={},
            unc="num",
            overflow=overflow,
        )
    rax.set_ylabel("$\dfrac{data}{Pred.}$")
    rax.set_ylim(0.5, 1.5)
    rax.set_xlim(xlimits)
    
        ## set axes labels and titles
    rax.set_xlabel(xlabel)

    return ax, rax


def plot_sys_variations(ax, rax, nom, up, dw, xlabel="", ylabel="", xlimits=None, ylimits=None, **opts):

    nom_mc = nom[mc_samples].integrate("process")
    nom_data = nom[data_samples].integrate("process")
    up_mc = up[mc_samples].integrate("process")
    dw_mc = dw[mc_samples].integrate("process")

    nom_opts = opts.get("Nominal", hstyles["Nominal"])
    up_opts = opts.get("Up", hstyles["Up"])
    dw_opts = opts.get("Down", hstyles["Down"])

    ## plot yields
        ## plot MC
    ax = hep.plot.histplot(nom_mc.values()[()], nom_mc.dense_axes()[0].edges(), ax=ax, label=nom_opts["name"], linestyle="-", color=nom_opts["color"], histtype="step")
    ax = hep.plot.histplot(up_mc.values()[()], up_mc.dense_axes()[0].edges(), ax=ax, label=up_opts["name"], linestyle="-", color=up_opts["color"], histtype="step")
    ax = hep.plot.histplot(dw_mc.values()[()], dw_mc.dense_axes()[0].edges(), ax=ax, label=dw_opts["name"], linestyle="-", color=dw_opts["color"], histtype="step")
        ## plot data
    ax = hep.plot.histplot(nom_data.values()[()], nom_data.dense_axes()[0].edges(), ax=ax, label="data", marker=".", color="k", histtype="errorbar")

    ax.legend(loc="upper right", title=opts["legend_title"])
    #ax.autoscale(axis="x", tight=True)
    ax.autoscale()
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, None)
    ax.set_xlabel(None)
    ax.set_xlim(xlimits)

    logy = opts.get("logy", False)
    if logy:
        nom_data_yrange = get_ylimits(nom_data, x_range=xlimits)
        nom_mc_yrange = get_ylimits(nom_mc, x_range=xlimits)
        up_mc_yrange = get_ylimits(up_mc, x_range=xlimits)
        dw_mc_yrange = get_ylimits(dw_mc, x_range=xlimits)
        yrange_vals = np.vstack((nom_mc_yrange, up_mc_yrange, dw_mc_yrange, nom_data_yrange))
        ax.set_yscale("log")
        #set_trace()
        ax.set_ylim(np.maximum(np.min(yrange_vals[:, 0]), 1e-1), 1.5*10**(np.ceil(np.log10(np.max(yrange_vals[:, 1])))))
    
        ## plot ratios
    plot.plotratio(
        nom_data, nom_mc, error_opts={"linestyle" : "-", "color" : nom_opts["color"]},
        unc="num", clear=False, ax=rax, guide_opts={},
        #unc="num", clear=False, ax=rax, denom_fill_opts={}, guide_opts={},
    )
    plot.plotratio(
        nom_data, up_mc, error_opts={"linestyle" : "-", "color" : up_opts["color"]},
        unc="num", clear=False, ax=rax, guide_opts={},
        #unc="num", clear=False, ax=rax, denom_fill_opts={}, guide_opts={},
    )
    plot.plotratio(
        nom_data, dw_mc, error_opts={"linestyle" : "-", "color" : dw_opts["color"]},
        unc="num", clear=False, ax=rax, guide_opts={},
        #unc="num", clear=False, ax=rax, denom_fill_opts={}, guide_opts={},
    )
    #rax.set_ylabel("data/MC")
    rax.set_ylabel("$\dfrac{data}{Pred.}$")
    rax.set_ylim(0.5, 1.5)
    rax.set_xlim(xlimits)
    rax.set_xlabel(xlabel)
    
    return ax, rax



def plot_mc1d(ax, hdict, xlabel="", ylabel="Events", xlimits=None, ylimits=None, hist_styles=None, **mc_opts):

    #set_trace()
    mcorder = mc_opts.get("mcorder")
    stack = mc_opts.get("stack", True)
    normalize = mc_opts.get("normalize", False)
    err_opts = mc_opts.get("error_opts", stack_error_opts)
    overflow = mc_opts.get("overflow", "none")
    logy = mc_opts.get("logy", False)
    leg_ncols = mc_opts.get("leg_ncols", 1)

        ## plot MC and data
    plot.plot1d(hdict[mc_samples],
        overlay=hdict.axes()[0].name,
        ax=ax,
        clear=False,
        stack=stack,
        line_opts=None,
        fill_opts=stack_fill_opts,
        error_opts=err_opts,
        order=mcorder,
        density=normalize,
        overflow=overflow,
        overlay_overflow=overflow,
    )
    #ax.autoscale(axis="x", tight=True)
    if logy:
        ax.set_yscale("log")
    ax.autoscale()
    ylims = (0, None) if ylimits is None else ylimits
    ax.set_ylim(ylims)
    #ax.set_ylim(0, ax.get_ylim()[1]*1.3)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_xlim(xlimits)

        ## set legend and corresponding colors
    handles, labels = ax.get_legend_handles_labels()
    for idx, sample in enumerate(labels):
        if sample == "data" or sample == "Observed": continue
        facecolor, legname = plt_tools.get_styles(sample, hstyles) if hist_styles is None else plt_tools.get_styles(sample, hist_styles)
        handles[idx].set_facecolor(facecolor)
        labels[idx] = legname
    # call ax.legend() with the new values
    #labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    ax.legend(handles,labels, loc="upper right", ncol=leg_ncols)
    
    return ax

def plot_2d_norm(values, xlimits, ylimits, xlabel="", ylabel="", hdict=None, mask=None, ax=None, xaxis_name=None, yaxis_name=None, xbins=None, ybins=None, **opts):
    if ax is None:
        ax = plt.gca()

    if xbins is not None:
        xedges = xbins
    else:
        xaxis = hdict.axis(xaxis_name)
        xedges = xaxis.edges()
    if ybins is not None:
        yedges = ybins
    else:
        yaxis = hdict.axis(yaxis_name)
        yedges = yaxis.edges()

    if mask is not None:
        values = np.ma.masked_where(mask, values)

    cmap = opts.get("cmap", "viridis")
    pc = ax.pcolormesh(xedges, yedges, values.T, cmap=cmap)
    ax.add_collection(pc)

    cmap_label = opts.get("cmap_label", "")
    plt.colorbar(pc, ax=ax, label=cmap_label, pad=0.)

    ax.autoscale()
    ax.set_xlim(xlimits)
    ax.set_ylim(ylimits)
    xtitle = opts.get("xtitle", xlabel)
    ytitle = opts.get("ytitle", ylabel)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    
    return ax


def plot_1D(values, bins, xlimits=None, xlabel="", ylabel="Events", linestyle="-", color="k", density=False, weights=None, ax=None, label="", histtype="step", show_step=True, sumw2=None, **kwargs):
    if ax is None:
        ax = plt.gca()

    if histtype == "errorbar":
        #set_trace()
        if sumw2 is None: raise ValueError("Sumw2 must be specified when plotting errorbars")
        hep.plot.histplot(values, bins, w2=sumw2, density=density, ax=ax, label=label, marker=None, color=color, histtype=histtype, xerr=True)
        #hep.plot.histplot(values, bins, density=density, ax=ax, label=label, marker=None, color=color, histtype=histtype, xerr=True)
        if show_step:
            hep.plot.histplot(values, bins, density=density, ax=ax, label=None, linestyle=linestyle, color=color, histtype="step")
    else:
        hep.plot.histplot(values, bins, density=density, ax=ax, label=label, linestyle=linestyle, color=color, histtype=histtype)
    #ax.autoscale(axis="x", tight=True)
    ax.autoscale()
    ax.set_xlim(min(bins), max(bins)) if xlimits is None else ax.set_xlim(xlimits)
    xtitle = kwargs.get("xtitle", xlabel)
    ytitle = kwargs.get("ytitle", ylabel)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    #ax.set_ylim(0, None)

    return ax


def QCD_Est(sig_reg, iso_sb, btag_sb, double_sb, norm_type=None, shape_region=None, norm_region=None, sys="nosys"):
    if not norm_type:
        raise ValueError("Normalization type has to be specified for qcd estimation")
    if not shape_region:
        raise ValueError("Region to get qcd shape has to be specified for qcd estimation")

    #set_trace()
    if sig_reg.sparse_dim() > 1:
        if sys in systematics.ttJets_sys.values():
            tt_dict = sig_reg["TT", sys].integrate("sys") if "TT" in sorted(set([key[0] for key in sig_reg.values().keys()])) else sig_reg["ttJets*", sys].integrate("sys") # added if statement when making template files
            non_tt_mc_dict = sig_reg[nonTT_mc_mask, "nosys"].integrate("sys")
            mc_dict = iso_sb.copy()
            mc_dict.clear()
            mc_dict.add(tt_dict)
            mc_dict.add(non_tt_mc_dict)
        else:
            mc_dict = sig_reg[:, sys].integrate("sys")
            mc_dict = mc_dict[mc_samples]
        data_dict = sig_reg[:, "nosys"].integrate("sys")
        data_dict = data_dict[data_samples]
            ## subsitute values into new sig_reg hist
        sig_reg = iso_sb.copy()
        sig_reg.clear()
        sig_reg.add(mc_dict)
        sig_reg.add(data_dict)
                
    
    sig_dmp = data_minus_prompt(sig_reg)
    iso_dmp = data_minus_prompt(iso_sb)
    btag_dmp = data_minus_prompt(btag_sb)
    double_dmp = data_minus_prompt(double_sb)

    data_histo_dict = {
        "ISO" : iso_sb[data_samples].integrate("process"),
        "BTAG" : btag_sb[data_samples].integrate("process"),
        "DOUBLE" : double_sb[data_samples].integrate("process"),
        "SIG" : sig_reg[data_samples].integrate("process"),
    }
    prompt_histo_dict = {
        "ISO" : iso_sb[prompt_mc_mask].integrate("process"),
        "BTAG" : btag_sb[prompt_mc_mask].integrate("process"),
        "DOUBLE" : double_sb[prompt_mc_mask].integrate("process"),
        "SIG" : sig_reg[prompt_mc_mask].integrate("process"),
    }
    dmp_dict = {
        "ISO" : iso_dmp,
        "BTAG" : btag_dmp,
        "DOUBLE" : double_dmp,
        "SIG" : sig_dmp,
    }

    qcd_dict = {
        "ISO" : iso_sb[qcd_samples].integrate("process"),
        "BTAG" : btag_sb[qcd_samples].integrate("process"),
        "DOUBLE" : double_sb[qcd_samples].integrate("process"),
        "SIG" : sig_reg[qcd_samples].integrate("process"),
    }

        # get normalized qcd shape (bins < 0. not allowed)
    qcd_norm_sumw, qcd_norm_sumw2 = get_qcd_shape(dmp_dict[shape_region])

    # find normalization
    """
    sideband norm = (QCD MC in sig region)/(QCD MC in norm_region)*(data-prompt in norm_region)

    ABCD norm = (BTAG*ISO)/(DOUBLE) data-prompt yields
    """
    normalization = 0
    if norm_type == "ABCD":
        normalization = (np.sum(btag_dmp.values(overflow="all")[()])*np.sum(iso_dmp.values(overflow="all")[()]))/(np.sum(double_dmp.values(overflow="all")[()]))
    elif norm_type == "ABCD-ML":
        normalization, norm_err = get_ml_qcd_yield(data=data_histo_dict, prompt=prompt_histo_dict, qcd=dmp_dict)
    elif norm_type == "Sideband":
        normalization = np.sum(qcd_dict["SIG"].values(overflow="all")[()])/np.sum(qcd_dict[norm_region].values(overflow="all")[()])*np.sum(dmp_dict[norm_region].values(overflow="all")[()])

        ## rescale yields and errors
    qcd_est_sumw = qcd_norm_sumw*normalization
    qcd_est_sumw2 = qcd_norm_sumw2*(normalization**2)

    output_qcd = sig_reg.copy()
        # substitute qcd_est array into original hist
    for idx in range(len(qcd_est_sumw2)):
        output_qcd.values(overflow="all", sumw2=True)[("QCD",)][0][idx] = qcd_est_sumw[idx]
        output_qcd.values(overflow="all", sumw2=True)[("QCD",)][1][idx] = qcd_est_sumw2[idx]
    
    return output_qcd



#def BKG_Est(sig_reg, sb_reg, norm_type=None, sys="nosys", ignore_uncs=False, isForTemplates=False):
def BKG_Est(sig_reg, sb_reg, norm_type=None, sys="nosys", isForTemplates=False, uncs_percentage=None):
    #set_trace()
    allowed_norm_types = ["SigMC", "Sideband"]
    #if norm_type != "SigMC":
    if norm_type not in allowed_norm_types:
        raise ValueError("Only SigMC supported for normalization type for bkg estimation right now.")

    if sig_reg.sparse_dim() > 1:
        if sys in systematics.ttJets_sys.values():
            tt_dict = sig_reg["TT", sys].integrate("sys") if "TT" in sorted(set([key[0] for key in sig_reg.values().keys()])) else sig_reg["ttJets*", sys].integrate("sys") # added if statement when making template files
            non_tt_mc_dict = sig_reg[nonTT_mc_mask, "nosys"].integrate("sys")
            mc_dict = iso_sb.copy()
            mc_dict.clear()
            mc_dict.add(tt_dict)
            mc_dict.add(non_tt_mc_dict)
        else:
            mc_dict = sig_reg[:, sys].integrate("sys")
            mc_dict = mc_dict[mc_samples]
        data_dict = sig_reg[:, "nosys"].integrate("sys")
        data_dict = data_dict[data_samples]
            ## subsitute values into new sig_reg hist
        #sig_reg = sb_reg.copy()
        sig_reg = mc_dict.copy()
        sig_reg.clear()
        sig_reg.add(mc_dict)
        sig_reg.add(data_dict)

        # data - (st + ttbar)
    dmt_dict = {
        "SB" : data_minus_top(sb_reg, isForTemplates),
        "SIG" : data_minus_top(sig_reg, isForTemplates),
    }
        # (EWK + QCD) contributions
    bkg_dict = {
        "SB" : sb_reg[bkg_mask].integrate("process"),
        "SIG" : sig_reg[bkg_mask].integrate("process"),
    }

        # get normalized bkg shape (bins < 0. not allowed)
    bkg_norm_sumw, bkg_norm_sumw2 = get_qcd_shape(dmt_dict["SB"])

    # find normalization
    """
    SigMC norm: yield from signal region MC
    """
    norm_sumw, norm_sumw2 = 0, 0
    if norm_type == "SigMC":
        #set_trace()
        norm_sumw, norm_sumw2 = np.sum(bkg_dict["SIG"].values(overflow="all", sumw2=True)[()], axis=1)
    elif norm_type == "Sideband":
        #set_trace()
        norm_sumw, norm_sumw2 = np.sum(bkg_dict["SIG"].values(overflow="all", sumw2=True)[()], axis=1)/np.sum(bkg_dict["SB"].values(overflow="all", sumw2=True)[()], axis=1)*np.sum(dmt_dict["SB"].values(overflow="all", sumw2=True)[()], axis=1)

        # rescale yields and errors
    bkg_est_sumw = bkg_norm_sumw * norm_sumw
    bkg_est_sumw2 = (bkg_est_sumw * (uncs_percentage/100))**2 if uncs_percentage else bkg_norm_sumw*norm_sumw2
    #set_trace()

    output_bkg = sig_reg.copy()
        # substitute qcd_est array into original hist
    output_bkg.values(overflow="all", sumw2=True)[("EWQCD",)][0][:] = bkg_est_sumw
    output_bkg.values(overflow="all", sumw2=True)[("EWQCD",)][1][:] = bkg_est_sumw2

    return output_bkg



def data_minus_prompt(histo):
    data_hist = histo[data_samples].integrate("process")
    pmc_hist  = histo[prompt_mc_mask].integrate("process")

    neg_pmc_hist = pmc_hist.copy()
    data_minus_prompt = data_hist.copy()
    neg_pmc_hist.scale(-1.)
    data_minus_prompt.add(neg_pmc_hist)

    return data_minus_prompt

def data_minus_top(histo, isForTemplates):
    #set_trace()
    data_hist = histo[data_samples].integrate("process")
    top_mc_hist  = histo[template_top_samples].integrate("process") if isForTemplates else histo[top_samples].integrate("process")

    neg_top_hist = top_mc_hist.copy()
    data_minus_top = data_hist.copy()
    neg_top_hist.scale(-1.)
    data_minus_top.add(neg_top_hist)

    return data_minus_top


def get_qcd_shape(dmp_hist, overflow="all"):
    sumw, sumw2 = dmp_hist.values(overflow=overflow, sumw2=True)[()]
        ## normalize shape
    sumw[sumw < 0] = 0
    norm_sumw = sumw/np.sum(sumw)
        ## "normalize" uncertainties
    norm_sumw2 = sumw2/(np.sum(sumw)**2)
    norm_sumw2[~np.isfinite(norm_sumw2)] = 0

    return norm_sumw, norm_sumw2



def Divide(h1, h2):
    set_trace()
    val1, val2 = h1.values(), h2.values()
    var1, var2 = h1.variances(), h2.variances()
    ratio_val = val1/val2
        # variance of 2 ratios estimated from Eq 4.12 assuming independent vars https://www.sagepub.com/sites/default/files/upm-binaries/6427_Chapter_4__Lee_(Analyzing)_I_PDF_6.pdf
    ratio_var = ratio_val**2 * ( (var1/val1**2) + (var2/val2**2) )

    hratio = h1.copy()
    hratio.values()[:], hratio.variances()[:] = ratio_val, ratio_var

    return hratio



def get_samples_yield_and_frac(histo, lumi, promptmc=False, overflow="all", sys="nosys"):
    """
    Get the yield and relative fraction for each sample of MC, get data yield and compare data/MC
    Returns: list of tuples containing sample name, yield, and fraction
    """
    if histo.sparse_dim() > 1:
        set_trace()
        mc_dict = histo[:, sys].integrate("sys")
        sig_dict = mc_dict[signal_samples]
        mc_dict = mc_dict[nonsignal_samples][mc_samples]
        data_dict = histo[:, "nosys"].integrate("sys")
        data_dict = data_dict[data_samples]

    else:
        sig_dict = histo[signal_samples]
        mc_dict = histo[nonsignal_samples][mc_samples]
        data_dict = histo[data_samples]

        # get signal values
    sig_yields_and_errs = sig_dict.integrate(sig_dict.dense_axes()[0].name).sum().values(overflow=overflow, sumw2=True)
    sig_proc_yields_list = [("".join(process), proc_yields[0], proc_yields[1]) for process, proc_yields in sig_yields_and_errs.items()]

        # get SM MC and data values
    yields_and_errs = mc_dict.integrate(mc_dict.dense_axes()[0].name).sum().values(overflow=overflow, sumw2=True)
    yields_and_errs.update(data_dict.integrate(data_dict.dense_axes()[0].name).sum().values(overflow=overflow, sumw2=True))
    proc_yields_list = [("".join(process), proc_yields[0], proc_yields[1]) for process, proc_yields in yields_and_errs.items()]
    mc_yield = sum([process[1] for process in proc_yields_list if not (process[0] == "data")])
    data_yield = sum([process[1] for process in proc_yields_list if (process[0] == "data")])
    mc_err = np.sqrt(sum([process[2] for process in proc_yields_list if not (process[0] == "data")]))
    data_err = np.sqrt(sum([process[2] for process in proc_yields_list if (process[0] == "data")]))

    rows = [("Lumi: %s fb^-1" % format(lumi, ".1f"), "Sample", "Yield", "Error", "Frac")]
    rows += [("", process, format(proc_yield, ".1f"), format(np.sqrt(proc_err), ".1f"), format((proc_yield/mc_yield)*100, ".1f")) for process, proc_yield, proc_err in proc_yields_list if not process == "data"]
    rows += [("", "SIM", format(mc_yield, ".1f"), format(mc_err, ".1f"), "100.0")]
    rows += [("", "", "", "", "")]
    rows += [("", "data", format(data_yield, ".1f"), format(data_err, ".1f"), "")]
    rows += [("", "", "", "", "")]
    rows += [("", "data/SIM", "", "", "0." if mc_yield == 0. else format(data_yield/mc_yield, ".3f"))]
    if promptmc:
        prompt_mc_yield = sum([process[1] for process in proc_yields_list if not ((process[0] == "data") or (process[0] == "QCD") or (process[0] == "EWK"))])
        prompt_mc_err = np.sqrt(sum([process[2] for process in proc_yields_list if not ((process[0] == "data") or (process[0] == "QCD") or (process[0] == "EWK"))]))
        rows += [("", "", "", "", "")]
        rows += [("", "Prompt MC", format(prompt_mc_yield, ".1f"), format(prompt_mc_err, ".1f"),  "")]
        rows += [("", "data-Prompt MC", format(data_yield-prompt_mc_yield, ".1f"), format(np.sqrt(data_err**2+prompt_mc_err**2), ".1f"), "")]

    rows += [("", "", "", "", "")]
    rows += [("", process, format(proc_yield, ".1f"), format(np.sqrt(proc_err), ".1f"), "") for process, proc_yield, proc_err in sig_proc_yields_list]

    yields_dict = {process:(round(proc_yield, 2), round(np.sqrt(proc_err), 2)) for process, proc_yield, proc_err in proc_yields_list}
    yields_dict.update({"SIM": (round(mc_yield, 2), round(mc_err, 2))})
    yields_dict.update({"data/SIM": 0. if mc_yield == 0. else round(data_yield/mc_yield, 3)})
    yields_dict.update({process:(round(proc_yield, 2), round(np.sqrt(proc_err), 2)) for process, proc_yield, proc_err in sig_proc_yields_list})

    return rows, yields_dict


def get_sys_variations_yield_and_frac(nosys, up, down, up_sys, dw_sys, overflow="all"):
    """
    Get the data/MC fraction and sys MC yield/nosys MC yield for nosys and up/down variations
    Returns: list of tuples containing systematics, data/MC fraction, and sys MC/nosys MC fraction
    """
        # get data yield from nosys (all same)
    data = np.sum(nosys[data_samples].integrate("process").values(overflow=overflow)[()])

        # get total MC yield
    nosys_mc = np.sum(nosys[mc_samples].integrate("process").values(overflow=overflow)[()]) # nosys
    up_mc = np.sum(up[mc_samples].integrate("process").values(overflow=overflow)[()]) # up variation
    down_mc = np.sum(down[mc_samples].integrate("process").values(overflow=overflow)[()]) # down variation

    rows = [("Systematic", "data/MC", "sys MC/nosys MC")]
    rows += [("nosys", format(data/nosys_mc, ".3f"), format(nosys_mc/nosys_mc, ".3f"))]
    rows += [("%s Up" % up_sys, format(data/up_mc, ".3f"), format(up_mc/nosys_mc, ".3f"))]
    rows += [("%s Down" % dw_sys, format(data/down_mc, ".3f"), format(down_mc/nosys_mc, ".3f"))]

    return rows



def get_ylimits(histo, x_range, mask_zero=True):
    ## get bin numbers corresponding to lower and upper values of x_range
    xmin_bin = np.where(histo.dense_axes()[0].edges() >= x_range[0])[0][0]
    xmax_bin = np.where(histo.dense_axes()[0].edges() >= x_range[1])[0][0]

    vals_array = histo.values()[()][xmin_bin:xmax_bin]

    min_yval = np.min(vals_array[vals_array > 0]) if mask_zero else np.min(vals_array) 
    max_yval = np.max(vals_array)

    return np.array([min_yval, max_yval])


def linearize_hist(histo, overflow=False, no_transpose=False, debug=False):
    from coffea import hist

    if histo.dense_dim() != 2:
        raise ValueError("Hist must be 2D in order to linearize!")
    if histo.sparse_dim() > 2:
        raise ValueError("Hist can have at most 2 sparse axes!")

    if debug: set_trace()
    xaxis = histo.dense_axes()[0]
    nbinsx = len(xaxis.edges())-1
    yaxis = histo.dense_axes()[1]
    nbinsy = len(yaxis.edges())-1
    nbins = nbinsx*nbinsy
    if debug: set_trace()
    output_hist = hist.Hist(
        "Events",
        *histo.sparse_axes(),
        hist.Bin("x_y","x_y", nbins, 0., nbins)
        #hist.Bin("%s_%s" % (xaxis.name, yaxis.name),"%s_%s" % (xaxis.name, yaxis.name), nbins, 0., nbins)
    )
    output_hist.axis("x_y")._bin_names = np.tile(histo.axis(xaxis).edges(), histo.axis(yaxis).edges().size)

        ## initialize hist to have empty bins for each key in histo.values().keys()
    for key in histo.values().keys():
        if len(key) == 0:
            output_hist.fill(x_y=np.zeros(0), weight=np.zeros(0))
        elif len(key) == 1:
            output_hist.fill(process=key[0], x_y=np.zeros(0), weight=np.zeros(0))
            #output_hist.fill(process=key[0], mtt_ctstar_abs=np.zeros(0), weight=np.zeros(0))
        else:
            output_hist.fill(process=key[0], sys=key[1], x_y=np.zeros(0), weight=np.zeros(0))
            #output_hist.fill(process=key[0], sys=key[1], mtt_ctstar_abs=np.zeros(0), weight=np.zeros(0))

        sumw_2D, sumw2_2D = histo.values(sumw2=True, overflow="all")[key] if overflow else histo.values(sumw2=True)[key]
        sumw_array = sumw_2D if no_transpose else sumw_2D.T
        sumw2_array = sumw2_2D if no_transpose else sumw2_2D.T

        if debug: set_trace()
        if overflow:
                ## combine over and under flow bins from ctstar dim with bins next to them
            sumw_comb_ct_overflow = np.zeros((sumw_array.shape[0]-2, sumw_array.shape[1]))
            sumw_comb_ct_overflow[0] = np.add(sumw_array[0], sumw_array[1])
            sumw_comb_ct_overflow[-1] = np.add(sumw_array[-2], sumw_array[-1])
            sumw_comb_ct_overflow[1:-1] = sumw_array[2:-2]
    
            sumw2_comb_ct_overflow = np.zeros((sumw2_array.shape[0]-2, sumw2_array.shape[1]))
            sumw2_comb_ct_overflow[0] = np.add(sumw2_array[0], sumw2_array[1])
            sumw2_comb_ct_overflow[-1] = np.add(sumw2_array[-2], sumw2_array[-1])
            sumw2_comb_ct_overflow[1:-1] = sumw2_array[2:-2]
    
                ## combine over and under flow bins from mtt dim with bins next to them
            sumw_comb_mtt_overflow = np.zeros((nbinsx, nbinsy))
            sumw_comb_mtt_overflow[0] = np.add(sumw_comb_ct_overflow.T[0], sumw_comb_ct_overflow.T[1])
            sumw_comb_mtt_overflow[-1] = np.add(sumw_comb_ct_overflow.T[-2], sumw_comb_ct_overflow.T[-1])
            sumw_comb_mtt_overflow[1:-1] = sumw_comb_ct_overflow.T[2:-2]
    
            sumw2_comb_mtt_overflow = np.zeros((nbinsx, nbinsy))
            sumw2_comb_mtt_overflow[0] = np.add(sumw2_comb_ct_overflow.T[0], sumw2_comb_ct_overflow.T[1])
            sumw2_comb_mtt_overflow[-1] = np.add(sumw2_comb_ct_overflow.T[-2], sumw2_comb_ct_overflow.T[-1])
            sumw2_comb_mtt_overflow[1:-1] = sumw2_comb_ct_overflow.T[2:-2]
    
                ## fill bins
            for ybin in range(nbinsy):
                output_hist.values(sumw2=True)[key][0][nbinsx*ybin:nbinsx*(ybin+1)] = sumw_comb_mtt_overflow.T[ybin]
                output_hist.values(sumw2=True)[key][1][nbinsx*ybin:nbinsx*(ybin+1)] = sumw2_comb_mtt_overflow.T[ybin]

        else:
                ## fill bins
            for ybin in range(nbinsy):
                output_hist.values(sumw2=True)[key][0][nbinsx*ybin:nbinsx*(ybin+1)] = sumw_array[ybin]
                output_hist.values(sumw2=True)[key][1][nbinsx*ybin:nbinsx*(ybin+1)] = sumw2_array[ybin]

    return output_hist



def root_converters_dict_to_hist(dict, vars=[], sparse_axes_list=[], dense_axes_list=[], overflow=False, transpose_da=False):
    """
    Inputs:
        dict : dictionary made from coffea "convert_histo_root_file"
        vars : list of variables to make histograms from
        sparse_axes_list : list of dictionaries with info to create sparse axes

    Output:
        dictionary of {var: coffea histogram objects}, each with 1 dense axis
    """
    from coffea import hist

    #set_trace()
    output_histos = {}
    for var in vars:
        if not (var, "dense_lookup") in dict.keys():
            raise ValueError(f"{var} not found in dense_lookup_dict")
        sumw = dict[(var, "dense_lookup")][0]
        sumw2 = dict[(f"{var}_error", "dense_lookup")][0]
        #edges = dict[(var, "dense_lookup")][1][0]
        edges = dict[(var, "dense_lookup")][1]

        sparse_axes = [hist.Cat(sparse_axis["name"], sparse_axis["label"]) for sparse_axis in sparse_axes_list]
        dense_axes = [hist.Bin(dense_axis["name"], dense_axis["name"], edges[dense_axis["idx"]]) for dense_axis in dense_axes_list]
        output_hist = hist.Hist(
            "Events",
            *sparse_axes,
            *dense_axes,
        )

        fill_dict = {sparse_axis["name"] : sparse_axis["fill"] for sparse_axis in sparse_axes_list}
        fill_dict.update({dense_axis["name"] : np.zeros(0) for dense_axis in dense_axes_list})
        fill_dict.update({"weight" : np.zeros(0)})
            ## initialize output_hist to have empty bins
        output_hist.fill(**fill_dict)

        #set_trace()
            ## fill bins
        for cat in output_hist.values().keys():
            if len(edges) == 1:
                output_hist.values(sumw2=True)[cat][0][:] = sumw
                output_hist.values(sumw2=True)[cat][1][:] = sumw2
            elif len(edges) == 2:
                for xbin in range(len(edges[0])-1):
                    output_hist.values(sumw2=True)[cat][0][xbin] = sumw.T[xbin] if transpose_da else sumw[xbin]
                    output_hist.values(sumw2=True)[cat][1][xbin] = sumw2.T[xbin] if transpose_da else sumw2[xbin]
            else: raise ValueError("Only 1D and 2D hists supported")

        output_histos[var] = output_hist

    return output_histos


def np_array_TO_hist(sumw, sumw2, hist_template, cat=None, overflow="none"):
    out_histo = hist_template.copy()
    out_histo.clear()

    fill_dict = {out_histo.sparse_axes()[0].name: cat} if len(out_histo.sparse_axes()) == 1 else {}
    fill_dict.update({out_histo.dense_axes()[idx].name: np.zeros(0) for idx in range(len(out_histo.dense_axes()))})
    fill_dict.update({"weight" : np.zeros(0)})
    out_histo.fill(**fill_dict)

        ## fill bins
    if len(out_histo.dense_axes()) == 1:
        if cat is None:
            for xbin in range(sumw.size):
                out_histo.values(overflow=overflow, sumw2=True)[()][0][xbin] = sumw[xbin]
                out_histo.values(overflow=overflow, sumw2=True)[()][1][xbin] = sumw2[xbin]
        else:
            for xbin in range(sumw.size):
                out_histo.values(overflow=overflow, sumw2=True)[(cat,)][0][xbin] = sumw[xbin]
                out_histo.values(overflow=overflow, sumw2=True)[(cat,)][1][xbin] = sumw2[xbin]
    elif len(out_histo.dense_axes()) == 2:
        if cat is None:
            for xbin in range(sumw.shape[0]):
                out_histo.values(overflow=overflow, sumw2=True)[()][0][xbin] = sumw[xbin]
                out_histo.values(overflow=overflow, sumw2=True)[()][1][xbin] = sumw2[xbin]
        else:
            for xbin in range(sumw.shape[0]):
                out_histo.values(overflow=overflow, sumw2=True)[(cat,)][0][xbin] = sumw[xbin]
                out_histo.values(overflow=overflow, sumw2=True)[(cat,)][1][xbin] = sumw2[xbin]
    else: raise ValueError("Only 1D and 2D hists supported")

    return out_histo 


def get_ratio_arrays(num_vals, denom_vals, input_bins):
    """
    Function that computes the ratio of two numpy arrays and returns the correctly formatted array to be used with step histogram
    """
    first_valid_bin, last_valid_bin = np.where(~np.isnan(num_vals/denom_vals))[0][0], np.where(~np.isnan(num_vals/denom_vals))[0][-1]+1
    ratio_masked_vals = np.ma.masked_where(np.isnan( (num_vals/denom_vals) [first_valid_bin:last_valid_bin]), (num_vals/denom_vals)[first_valid_bin:last_valid_bin])

    ratio_vals = np.r_[ratio_masked_vals, ratio_masked_vals[-1]]
    ratio_bins = input_bins[first_valid_bin:last_valid_bin+1]

    return ratio_vals, ratio_bins


## get uncertainty on ratio between A and B
def get_ratio_and_uncertainty(a_val, a_err, b_val, b_err):
        # c = a/b
    c_val = a_val/b_val
    c_err = np.sqrt( (a_err/a_val)**2 + (b_err/b_val)**2 )*c_val
    return c_val, c_err



def smoothing_mttbins(nosys, systematic, mtt_centers, nbinsx, nbinsy, debug=False, ratio=False, **fit_opts):
        # smoothing
    import statsmodels.nonparametric.smoothers_lowess as sm
    LOWESS = sm.lowess


    frac_val = fit_opts.get("frac", 0.5) # fraction of the data used when estimating each y-value
    it_val = fit_opts.get("it", 3)

    if debug: set_trace()
        # np array of original bin values
    nom_vals = nosys if isinstance(nosys, np.ndarray) else nosys.values()[()]
    systematic_vals = systematic if isinstance(systematic, np.ndarray) else systematic.values()[()]
    total_array = np.zeros(nbinsx*nbinsy)

        # loop over each bin of cos theta
    #set_trace()
    for ybin in range(nbinsy):
        yin = (systematic_vals[ybin*nbinsx:(ybin+1)*nbinsx] - nom_vals[ybin*nbinsx:(ybin+1)*nbinsx])/nom_vals[ybin*nbinsx:(ybin+1)*nbinsx] # relative deviation from nosys
        total_array[ybin*nbinsx:(ybin+1)*nbinsx] = np.nan_to_num(LOWESS(yin, mtt_centers, frac=frac_val, it=it_val, return_sorted=False)) # fill nan values with 0

    if debug: set_trace()
        # substitute smoothed array into copy of original hist
    #set_trace()
    if isinstance(systematic, np.ndarray):
        if ratio: # save sys/nominal
            smoothed_histo = 1+total_array
        else:
            smoothed_histo = (1+total_array)*nom_vals
    else:
        smoothed_histo = systematic.copy()
        for idx in range(nbinsx*nbinsy):
            if ratio: # save sys/nominal
                smoothed_histo.values()[()][idx] = 1+total_array[idx]
            else:
                smoothed_histo.values()[()][idx] = (1+total_array[idx])*nom_vals[idx]

    if debug: set_trace()
    return smoothed_histo


def new_smoothing(nosys, systematic, mtt_bins, nbinsx, nbinsy, debug=False, ratio=False, systematic_sumw2=None, p_end=0.01, cutoff=3.0, num_pt_comp=5):
        # smoothing
    import smoothing_compiled.pysmoother as pysmoother
    from coffea import hist
    if not isinstance(systematic, hist.Hist):
        if not isinstance(systematic_sumw2, np.ndarray):
            raise ValueError("Systematic must be coffea hist object or input sumw2 values for systematic variation")

    if debug: set_trace()
        # np array of original bin values
    nom_vals = nosys if isinstance(nosys, np.ndarray) else nosys.values()[()]
    if systematic_sumw2:
        sys_sumw = systematic.values()[()]
        sys_sumw = systematic_sumw2
    else:
        sys_sumw, sys_sumw2 = systematic.values(sumw2=True)[()]

        # get input arrays and values for py_smoothhist
    y_sumw = sys_sumw/nom_vals
    y_err = np.sqrt(sys_sumw2)/nom_vals
    total_array = np.zeros(nbinsx*nbinsy)

        # loop over each bin of cos theta
    #set_trace()
    for ybin in range(nbinsy):
        filtered_y_sumw = calculate_zscore(y_sumw[ybin*nbinsx:(ybin+1)*nbinsx], cutoff=cutoff, num_pt_comp=num_pt_comp)
        #set_trace()
        pysmoother.py_smoothhist(mtt_bins, np.nan_to_num(filtered_y_sumw), np.nan_to_num(y_err[ybin*nbinsx:(ybin+1)*nbinsx]), p_end, total_array[ybin*nbinsx:(ybin+1)*nbinsx])
        ##pysmoother.py_smoothhist(mtt_bins, np.nan_to_num(y_sumw[ybin*nbinsx:(ybin+1)*nbinsx]), np.nan_to_num(y_err[ybin*nbinsx:(ybin+1)*nbinsx]), p_end, total_array[ybin*nbinsx:(ybin+1)*nbinsx])

    if debug: set_trace()
        # substitute smoothed array into copy of original hist
    #set_trace()
    if isinstance(systematic, np.ndarray):
        if ratio: # save sys/nominal
            smoothed_histo = total_array
        else:
            smoothed_histo = total_array*nom_vals
    else:
        smoothed_histo = systematic.copy()
        for idx in range(nbinsx*nbinsy):
            if ratio: # save sys/nominal
                smoothed_histo.values()[()][idx] = total_array[idx]
            else:
                smoothed_histo.values()[()][idx] = total_array[idx]*nom_vals[idx]

    if debug: set_trace()
    return smoothed_histo




    # function that flattens systematic variations
def flatten(nosys, systematic, ratio=False):
        # np array of original bin values
    nom_vals = nosys.values()[()]
    systematic_vals = systematic.values()[()]

    flat_val= np.sum(systematic_vals)/np.sum(nom_vals)

        # substitute flattened array into copy of original hist
    flattened_histo = systematic.copy()
    for idx in range(nom_vals.size):
        if ratio:
            flattened_histo.values()[()][idx] = flat_val
        else:
            flattened_histo.values()[()][idx] = flat_val*nom_vals[idx]

        # save flat_val in overflow (in terms of relative deviation from nominal)
    flattened_histo.values(overflow="over")[()][-1] = flat_val - 1.0

    return flattened_histo



def calculate_zscore(vals, cutoff=3.0, num_pt_comp=5):
    #set_trace()
    #num_pt_comp = 5 # number of points on either side of point to be used for finding mean, stdev
    ##num_pt_comp = 3 # number of points on either side of point to be used for finding mean, stdev
    #cutoff = 3.0
    #cutoff = 4.0
    #cutoff = 2.0
    z_scores = -1*np.ones(vals.size)
    filtered_vals = vals.copy()

    for i in range(vals.size):
        sel_val = vals[i] # get selected value to excluded from z-score calculation
        if i < num_pt_comp: # indices on lower edge of values
            if i == 0: continue
            #if i == 0: # for first bin, only use half of specified num_pt_comp points instead of skipping
            #    bracket = vals[:int(np.ceil(num_pt_comp/2))+1]
            else:
                bracket = vals[:num_pt_comp+1]
        elif i >= vals.size - num_pt_comp: # indices on upper edge of values
            if i == vals.size - 1: # for last bin only use half of specified num_pt_comp points
                bracket = vals[-int(np.ceil(num_pt_comp/2))-1:]
            else:
                bracket = vals[-num_pt_comp-1:]
        else: # indices in middle
            bracket = vals[i-num_pt_comp:i+num_pt_comp+1]

        excluded = bracket[bracket != sel_val] # remove value corresponding to current index
        mean, std = excluded.mean(), excluded.std(ddof=1) #uses unbiased estimator for stdev
        z_scores[i] = np.nan if abs((sel_val - mean)/std) == np.inf else abs((sel_val - mean)/std)

        if z_scores[i] > cutoff:
            filtered_vals[i] = mean
    #set_trace()

    return filtered_vals

