"""
This script is a 2D extension of the coffea export1d function 
https://github.com/CoffeaTeam/coffea/blob/master/coffea/hist/export.py

Created 20 October 2021
"""

import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from uproot3_methods.classes.TH1 import Methods as TH1Methods
    from uproot3_methods.classes.TH2 import Methods as TH2Methods


class TH1(TH1Methods, list):
    pass

class TAxis(object):
    def __init__(self, fNbins, fXmin, fXmax):
        self._fNbins = fNbins
        self._fXmin = fXmin
        self._fXmax = fXmax

class TH2(TH2Methods, list):
    pass


def export1d(hist):
    """Export a 1-dimensional `Hist` object to uproot
    This allows one to write a coffea histogram into a ROOT file, via uproot.
    Parameters
    ----------
        hist : Hist
            A 1-dimensional histogram object
    Returns
    -------
        out
            A ``uproot3_methods.classes.TH1`` object
    Examples
    --------
    Creating a coffea histogram, filling, and writing to a file::

        import uproot3, numpy
        from coffea import hist
        h1D = hist.Hist("Events", hist.Bin("xaxis", "some x variable", 10, 0, 1))
        h1D.fill(xaxis=numpy.random.default_rng().uniform(low=0., high=1., size=1000))
        fout = uproot3.create("output.root")
        fout["1D"] = export1d(h1D)
        fout.close()
    """
    if hist.dense_dim() != 1:
        raise ValueError("export1d() can only support one dense dimension")
    if hist.sparse_dim() != 0:
        raise ValueError("export1d() expects zero sparse dimensions")

    axis = hist.axes()[0]
    sumw, sumw2 = hist.values(sumw2=True, overflow="all")[()]
    edges = axis.edges(overflow="none")

    out = TH1.__new__(TH1)
    out._fXaxis = TAxis(len(edges) - 1, edges[0], edges[-1])
    out._fXaxis._fName = axis.name
    out._fXaxis._fTitle = axis.label
    if not axis._uniform:
        out._fXaxis._fXbins = edges.astype(">f8")

    centers = (edges[:-1] + edges[1:]) / 2.0
    out._fEntries = out._fTsumw = out._fTsumw2 = sumw[1:-1].sum()
    out._fTsumwx = (sumw[1:-1] * centers).sum()
    out._fTsumwx2 = (sumw[1:-1] * centers ** 2).sum()

    out._fName = "histogram"
    out._fTitle = hist.label

    out._classname = b"TH1D"
    out.extend(sumw.astype(">f8"))
    out._fSumw2 = sumw2.astype(">f8")

    return out


def export2d(hist):
    """Export a 2-dimensional `Hist` object to uproot
    This allows one to write a coffea histogram into a ROOT file, via uproot.
    Parameters
    ----------
        hist : Hist
            A 2-dimensional histogram object
    Returns
    -------
        out
            A ``uproot3_methods.classes.TH2`` object
    Examples
    --------
    Creating a coffea histogram, filling, and writing to a file::
        import uproot3, numpy
        from coffea import hist
        h2D = hist.Hist("Events", hist.Bin("xaxis", "some x variable", 10, 0, 1), hist.Bin("yaxis", "some y variable", 20, 0, 10))
        h2D.fill(xaxis=numpy.random.default_rng().uniform(low=0., high=1., size=1000), yaxis=numpy.random.default_rng().uniform(low=0., high=10., size=1000))
        fout = uproot3.create("output.root")
        fout["2D"] = export2d(h2D)
        fout.close()
    """
    if hist.dense_dim() != 2:
        raise ValueError("export2d() can only support two dense dimension")
    if hist.sparse_dim() != 0:
        raise ValueError("export2d() expects zero sparse dimensions")

    xaxis, yaxis = hist.axes()[0], hist.axes()[1]
    sumw, sumw2 = hist.values(sumw2=True, overflow="all")[()]
    xedges, yedges = xaxis.edges(overflow="none"), yaxis.edges(overflow="none")

    out = TH2.__new__(TH2)
    out._fEntries = out._fTsumw = out._fTsumw2 = sumw.sum()
    # xaxis
    out._fXaxis = TAxis(len(xedges) - 1, xedges[0], xedges[-1])
    out._fXaxis._fName = xaxis.name
    out._fXaxis._fTitle = xaxis.label
    if not xaxis._uniform:
        out._fXaxis._fXbins = xedges.astype(">f8")
    # yaxis
    out._fYaxis = TAxis(len(yedges) - 1, yedges[0], yedges[-1])
    out._fYaxis._fName = yaxis.name
    out._fYaxis._fTitle = yaxis.label
    if not yaxis._uniform:
        out._fYaxis._fXbins = yedges.astype(">f8")

    # xaxis
    xcenters = (xedges[:-1] + xedges[1:]) / 2.
    out._fTsumwx = xcenters.dot(sumw[1:-1, 1:-1].sum(1))
    out._fTsumwx2 = (xcenters**2).dot(sumw[1:-1, 1:-1].sum(1))

    # yaxis
    ycenters = (yedges[:-1] + yedges[1:]) / 2.
    out._fTsumwy = ycenters.dot(sumw[1:-1, 1:-1].sum(0))
    out._fTsumwy2 = (ycenters**2).dot(sumw[1:-1, 1:-1].sum(0))

    out._fName = "histogram"
    out._fTitle = hist.label

    out._classname = b"TH2D"
    out.extend((sumw.T).astype(">f8"))
    out._fSumw2 = (sumw2.T).astype(">f8")

    return out
