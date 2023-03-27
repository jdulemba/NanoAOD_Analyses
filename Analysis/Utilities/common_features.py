#!/usr/bin/env python

"""
This script is meant to keep track of various things that are common to all histograms
or plotting features, such as bin/axis labels.
"""

jet_mults = {
    "3Jets" : "3 jets",
    "4PJets" : "4+ jets",
    "4Jets" : "4 jets",
    "5PJets" : "5+ jets",
}

objtypes = {
    "Jets" : "jets",
    "Muon" : "$\\mu$",
    "Electron" : "$e$",
    "Lepton" : "$\\ell$",
}

year_labels = {
    "2016APV" : "2016preVFP",
    "2016" : "2016postVFP",
    "2017" : "2017",
    "2018" : "2018",
    "Total" : "Run 2",
}

channel_labels = {
    "Muon_3Jets" : "$\\mu$/3 jets",
    "Muon_4PJets" : "$\\mu$/4+ jets",
    "Electron_3Jets" : "$e$/3 jets",
    "Electron_4PJets" : "$e$/4+ jets",
    "Lepton_3Jets" : "$\\ell$/3 jets",
    "Lepton_4PJets" : "$\\ell$/4+ jets",
    "Lepton_4Jets" : "$\\ell$/4 jets",
    "Lepton_5PJets" : "$\\ell$/5+ jets",
}

variable_names_to_labels = {
    "mtt" : "$m_{t\\bar{t}}$ [GeV]",
    "tlep_ctstar_abs" : "|$\cos(\\theta^{*}_{t_{\ell}})$|",
    "tlep_ctstar" : "$\cos(\\theta^{*}_{t_{\ell}})$",
    "mthad" : "$m_{t_{h}}$ [GeV]",
    "mWHad" : "$m_{W_{h}}$ [GeV]",
    "mWLep" : "$m_{W_{\ell}}$ [GeV]",
    "pt_thad" : "$p_{T}$($t_{h}$) [GeV]",
    "pt_tlep" : "$p_{T}$($t_{\ell}$) [GeV]",
    "pt_tt" : "$p_{T}$($t\\bar{t}$) [GeV]",
    "eta_thad" : "$\\eta$($t_{h}$)",
    "eta_tlep" : "$\\eta$($t_{\ell}$)",
    "eta_tt" : "$\\eta$($t\\bar{t}$)",
    "phi_tt" : "$\\phi$($t\\bar{t}$)",
    "full_disc" : "$\\lambda_{j}$",
    "mass_disc" : "-log($\\lambda_{Mj}$)",
    "ns_disc" : "-log($\\lambda_{NSj}$)",
    "ns_dist" : "$D_{\\nu, min}$ [GeV]",
    "Jets_njets" : "$n_{jets}$",
    "Jets_pt" : "$p_{T}$(jets) [GeV]",
    "Jets_eta" : "$\\eta$(jets)",
    "Jets_phi" : "$\\phi$(jets)",
    "Jets_LeadJet_pt" : "$p_{T}$(leading jet) [GeV]",
    "Jets_LeadJet_eta" : "$\\eta$(leading jet)",
    "Jets_DeepCSV_bDisc" : "DeepCSV b Disc",
    "Jets_DeepJet_bDisc" : "DeepJet b Disc",
    "Lep_pt" : "$p_{T}$(%s) [GeV]",
    "Lep_eta" : "$\\eta$(%s)",
    "Lep_phi" : "$\\phi$(%s)",
    "Lep_iso" : "pfRelIso, %s",
    "MT" : "$M_{T}$ [GeV]",
    "MET_pt" : "$p_{T}$(MET) [GeV]",
    "MET_phi" : "$\phi$(MET)",
    "ctstar" : "$\cos(\\theta^{*}_{t})$",
    "ctstar_abs" : "|$\cos(\\theta^{*}_{t})$|",
    "nusolver_chi2" : "$\\chi_{\\nu}^{2}$",
    "thad_ctstar_abs" : "|$\cos(\\theta^{*}_{t_{h}})$|",
    "thad_ctstar" : "$\cos(\\theta^{*}_{t_{h}})$",
    "rapidity_top" : "$y_t$",
    "mtop" : "$m_{t}$ [GeV]",
    "pt_top" : "$p_{T}$($t$) [GeV]",
    "eta_top" : "$\\eta$($t$)",
    "phi_top" : "$\\phi$($t$)",
    "top_ctstar_abs" : "|$\cos(\\theta^{*}_{t})$|",
    "top_ctstar" : "$\cos(\\theta^{*}_{t})$",
    "mtbar" : "$m_{\\bar{t}}$ [GeV]",
    "pt_tbar" : "$p_{T}$($\\bar{t}$) [GeV]",
    "eta_tbar" : "$\\eta$($\\bar{t}$)",
    "phi_tbar" : "$\\phi$($\\bar{t}$)",
    "tbar_ctstar_abs" : "|$\cos(\\theta^{*}_{\\bar{t}})$|",
    "tbar_ctstar" : "$\cos(\\theta^{*}_{\\bar{t}})$",
}    
