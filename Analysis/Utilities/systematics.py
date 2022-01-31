"""
List of systematics
"""

ind_jec_sys = {
    year : {
        #"JES_AbsoluteFlavMap",
        #"JES_CorrelationGroupFlavor",
        #"JES_CorrelationGroupIntercalibration",
        #"JES_CorrelationGroupMPFInSitu",
        #"JES_CorrelationGroupUncorrelated",
        #"JES_CorrelationGroupbJES",
        #"JES_FlavorPhotonJet",
        #"JES_FlavorPureBottom",
        #"JES_FlavorPureCharm",
        #"JES_FlavorPureGluon",
        #"JES_FlavorPureQuark",
        #"JES_FlavorZJet",
        #"JES_PileUpEnvelope",
        #"JES_PileUpMuZero",
        #"JES_PileUpPtEC2",
        #"JES_PileUpPtHF",
        #"JES_RelativeJEREC2",
        #"JES_RelativeJERHF",
        #"JES_RelativePtEC2",
        #"JES_RelativePtHF",
        #"JES_RelativeStatHF",
        #"JES_SubTotalAbsolute",
        #"JES_SubTotalMC",
        #"JES_SubTotalPileUp",
        #"JES_SubTotalPt",
        #"JES_SubTotalRelative",
        #"JES_SubTotalScale",
        #"JES_TotalNoFlavor",
        #"JES_TotalNoFlavorNoTime",
        #"JES_TotalNoTime",
        "JES_AbsoluteMPFBias",
        "JES_AbsoluteScale",
        "JES_AbsoluteStat",
        "JES_FlavorQCD",
        "JES_Fragmentation",
        "JES_PileUpDataMC",
        "JES_PileUpPtBB",
        "JES_PileUpPtEC1",
        "JES_PileUpPtRef",
        "JES_RelativeBal",
        "JES_RelativeFSR",
        "JES_RelativeJEREC1",
        "JES_RelativePtBB",
        "JES_RelativePtEC1",
        "JES_RelativeSample",
        "JES_RelativeStatEC",
        "JES_RelativeStatFSR",
        "JES_SinglePionECAL",
        "JES_SinglePionHCAL",
        "JES_TimePtEta",
        "JES_Total",
    }
    for year in ["2016APV", "2016", "2017", "2018"]
}


regroup_jec_sys = {
    year : {
        "JES_Absolute",
        f"JES_Absolute_{year}",
        "JES_BBEC1",
        f"JES_BBEC1_{year}",
        "JES_FlavorQCD",
        "JES_RelativeBal",
        f"JES_RelativeSample_{year}",
    }
    for year in ["2016APV", "2016", "2017", "2018"]
}


ttJets_sys = {
    "ISR_UP" : "ISRUp",
    "ISR_DW" : "ISRDown",
    "FSR_UP" : "FSRUp",
    "FSR_DW" : "FSRDown",
    "FACTOR_UP" : "FACTORUp",
    "FACTOR_DW" : "FACTORDown",
    "RENORM_UP" : "RENORMUp",
    "RENORM_DW" : "RENORMDown",
    "RENORM_FACTOR_UP" : "RENORM_FACTOR_SAMEUp",
    "RENORM_FACTOR_DW" : "RENORM_FACTOR_SAMEDown",
    "RENORM_UP_FACTOR_DW" : "RENORM_FACTOR_DIFFUp",
    "RENORM_DW_FACTOR_UP" : "RENORM_FACTOR_DIFFDown",
    "MTOP_UP" : "mtopUP",
    "MTOP_DW" : "mtopDOWN",
    "HDAMP_UP" : "hdampUP",
    "HDAMP_DW" : "hdampDOWN",
    "UE_UP" : "ueUP",
    "UE_DW" : "ueDOWN",
    "fsr_UP"  : "fsrUP",
    "fsr_DW": "fsrDOWN",
    "isr_UP"  : "isrUP",
    "isr_DW": "isrDOWN",
    "MTOP3GeV_DW": "mtop1695",
    "MTOP3GeV_UP": "mtop1755",
    
}

# dict of name_I_want : name_in_code
    # systematics that only change event weights, for nonsignal datasets
reweight_sys_opts = {
    year : {
        "LEP_RECO_UP" : "Lep_RECOUp",
        "LEP_RECO_DW" : "Lep_RECODown",
        "LEP_TRIG_UP" : "Lep_TRIGUp",
        "LEP_TRIG_DW" : "Lep_TRIGDown",
        "BTAG_BC_UP" : "btag_bc_up",
        "BTAG_BC_UP_CORR" : "btag_bc_up_correlated",
        "BTAG_BC_UP_UNCORR" : "btag_bc_up_uncorrelated",
        "BTAG_BC_DW" : "btag_bc_down",
        "BTAG_BC_DW_CORR" : "btag_bc_down_correlated",
        "BTAG_BC_DW_UNCORR" : "btag_bc_down_uncorrelated",
        "BTAG_L_UP" : "btag_l_up",
        "BTAG_L_UP_CORR" : "btag_l_up_correlated",
        "BTAG_L_UP_UNCORR" : "btag_l_up_uncorrelated",
        "BTAG_L_DW" : "btag_l_down",
        "BTAG_L_DW_CORR" : "btag_l_down_correlated",
        "BTAG_L_DW_UNCORR" : "btag_l_down_uncorrelated",
        #"BTAG_BC_UP" : "btag_bc_UP",
        #"BTAG_BC_DW" : "btag_bc_DW",
        #"BTAG_L_UP" : "btag_l_UP",
        #"BTAG_L_DW" : "btag_l_DW",
        "PILEUP_UP" : "PileupUp",
        "PILEUP_DW" : "PileupDown",
        "PREFIRE_UP" : "PrefireUp",
        "PREFIRE_DW" : "PrefireDown",
        "ISR_UP" : "ISRUp",
        "ISR_DW" : "ISRDown",
        "FSR_UP" : "FSRUp",
        "FSR_DW" : "FSRDown",
        "FACTOR_UP" : "FACTORUp",
        "FACTOR_DW" : "FACTORDown",
        "RENORM_UP" : "RENORMUp",
        "RENORM_DW" : "RENORMDown",
        "RENORM_FACTOR_UP" : "RENORM_FACTOR_SAMEUp",
        "RENORM_FACTOR_DW" : "RENORM_FACTOR_SAMEDown",
        "RENORM_UP_FACTOR_DW" : "RENORM_FACTOR_DIFFUp",
        "RENORM_DW_FACTOR_UP" : "RENORM_FACTOR_DIFFDown",
        "deltaQCDdeltaEW" : "EWcorrUp",
        #"EW_CORR_UP" : "EWcorrUp",
        "YUKAWA_UP" : "YukawaUp",
        "YUKAWA_DW" : "YukawaDown",
    }
    for year in ["2016APV", "2016", "2017", "2018"]
}
reweight_sys_opts["2018"].pop("PREFIRE_UP")
reweight_sys_opts["2018"].pop("PREFIRE_DW")

    # systematics that change entire object selection
event_sys_opts = {
    year : {
        "JER_UP" : "JER_UP",
        "JER_DW" : "JER_DW",
        "MET_UP" : "MET_UP",
        "MET_DW" : "MET_DW",
        #"JES_UP" : "JES_UP",
        #"JES_DW" : "JES_DW",
    }
    for year in ["2016APV", "2016", "2017", "2018"]
}
    # add individual JES systematics
event_sys_opts["2016APV"].update({"_".join([sys, "UP"]): "_".join([sys, "UP"]) for sys in ind_jec_sys["2016APV"]}) # add up
event_sys_opts["2016APV"].update({"_".join([sys, "DW"]): "_".join([sys, "DW"]) for sys in ind_jec_sys["2016APV"]}) # add down
event_sys_opts["2016"].update({"_".join([sys, "UP"]): "_".join([sys, "UP"]) for sys in ind_jec_sys["2016"]}) # add up
event_sys_opts["2016"].update({"_".join([sys, "DW"]): "_".join([sys, "DW"]) for sys in ind_jec_sys["2016"]}) # add down
event_sys_opts["2017"].update({"_".join([sys, "UP"]): "_".join([sys, "UP"]) for sys in ind_jec_sys["2017"]}) # add up
event_sys_opts["2017"].update({"_".join([sys, "DW"]): "_".join([sys, "DW"]) for sys in ind_jec_sys["2017"]}) # add down
event_sys_opts["2018"].update({"_".join([sys, "UP"]): "_".join([sys, "UP"]) for sys in ind_jec_sys["2018"]}) # add up
event_sys_opts["2018"].update({"_".join([sys, "DW"]): "_".join([sys, "DW"]) for sys in ind_jec_sys["2018"]}) # add down

#    # add regrouped JES systematics
#event_sys_opts["2016APV"].update({"_".join([sys, "UP"]): "_".join([sys, "UP"]) for sys in regroup_jec_sys["2016APV"]}) # add up
#event_sys_opts["2016APV"].update({"_".join([sys, "DW"]): "_".join([sys, "DW"]) for sys in regroup_jec_sys["2016APV"]}) # add down
#event_sys_opts["2016"].update({"_".join([sys, "UP"]): "_".join([sys, "UP"]) for sys in regroup_jec_sys["2016"]}) # add up
#event_sys_opts["2016"].update({"_".join([sys, "DW"]): "_".join([sys, "DW"]) for sys in regroup_jec_sys["2016"]}) # add down
#event_sys_opts["2017"].update({"_".join([sys, "UP"]): "_".join([sys, "UP"]) for sys in regroup_jec_sys["2017"]}) # add up
#event_sys_opts["2017"].update({"_".join([sys, "DW"]): "_".join([sys, "DW"]) for sys in regroup_jec_sys["2017"]}) # add down
#event_sys_opts["2018"].update({"_".join([sys, "UP"]): "_".join([sys, "UP"]) for sys in regroup_jec_sys["2018"]}) # add up
#event_sys_opts["2018"].update({"_".join([sys, "DW"]): "_".join([sys, "DW"]) for sys in regroup_jec_sys["2018"]}) # add down


# dict of name_I_want : name_in_code
    # systematics available to be applied to signal (reweighting)
signal_reweight_opts = {
    year : {
        "LEP_RECO_UP" : "Lep_RECOUp",
        "LEP_RECO_DW" : "Lep_RECODown",
        "LEP_TRIG_UP" : "Lep_TRIGUp",
        "LEP_TRIG_DW" : "Lep_TRIGDown",
        "BTAG_BC_UP" : "btag_bc_up",
        "BTAG_BC_UP_CORR" : "btag_bc_up_correlated",
        "BTAG_BC_UP_UNCORR" : "btag_bc_up_uncorrelated",
        "BTAG_BC_DW" : "btag_bc_down",
        "BTAG_BC_DW_CORR" : "btag_bc_down_correlated",
        "BTAG_BC_DW_UNCORR" : "btag_bc_down_uncorrelated",
        "BTAG_L_UP" : "btag_l_up",
        "BTAG_L_UP_CORR" : "btag_l_up_correlated",
        "BTAG_L_UP_UNCORR" : "btag_l_up_uncorrelated",
        "BTAG_L_DW" : "btag_l_down",
        "BTAG_L_DW_CORR" : "btag_l_down_correlated",
        "BTAG_L_DW_UNCORR" : "btag_l_down_uncorrelated",
        #"BTAG_BC_UP" : "btag_bc_UP",
        #"BTAG_BC_DW" : "btag_bc_DW",
        #"BTAG_L_UP" : "btag_l_UP",
        #"BTAG_L_DW" : "btag_l_DW",
        "PILEUP_UP" : "PileupUp",
        "PILEUP_DW" : "PileupDown",
        "PREFIRE_UP" : "PrefireUp",
        "PREFIRE_DW" : "PrefireDown",
        #"ISR_UP" : "ISRUp",
        #"ISR_DW" : "ISRDown",
        #"FSR_UP" : "FSRUp",
        #"FSR_DW" : "FSRDown",
        "AH_FACTOR_UP" : "AH_FACTORUp",
        "AH_FACTOR_DW" : "AH_FACTORDown",
        "AH_RENORM_UP" : "AH_RENORMUp",
        "AH_RENORM_DW" : "AH_RENORMDown",
        "AH_RENORM_FACTOR_UP" : "AH_RENORM_FACTOR_SAMEUp",
        "AH_RENORM_FACTOR_DW" : "AH_RENORM_FACTOR_SAMEDown",
        "AH_RENORM_UP_FACTOR_DW" : "AH_RENORM_FACTOR_DIFFUp",
        "AH_RENORM_DW_FACTOR_UP" : "AH_RENORM_FACTOR_DIFFDown",
    }
    for year in ["2016APV", "2016", "2017", "2018"]
}
signal_reweight_opts["2018"].pop("PREFIRE_UP")
signal_reweight_opts["2018"].pop("PREFIRE_DW")




# dict of systematic in analyzer to plotting name
sys_to_name = {
    year : {
        #"nosys" : "nosys",
        "EWcorrUp" : "deltaQCDdeltaEW",
        #"EWcorrUp" : "EW_CORR_UP",
        "YukawaUp" : "YUKAWA_UP",
        "YukawaDown" : "YUKAWA_DW",
        "shapeDown" : "EWQCD_SHAPE_DW",
        "shapeUp"   : "EWQCD_SHAPE_UP",
        "Lep_RECOUp"   : "LEP_RECO_UP",
        "Lep_RECODown" : "LEP_RECO_DW",
        "Lep_TRIGUp"   : "LEP_TRIG_UP",
        "Lep_TRIGDown" : "LEP_TRIG_DW",
        #"btag_bc_UP" : "BTAG_BC_UP",
        "btag_bc_up" : "BTAG_BC_UP",
        "btag_bc_up_correlated" : "BTAG_BC_CORR_UP",
        "btag_bc_up_uncorrelated" : "BTAG_BC_UNCORR_UP",
        #"btag_bc_DW" : "BTAG_BC_DW",
        "btag_bc_down" : "BTAG_BC_DW",
        "btag_bc_down_correlated" : "BTAG_BC_CORR_DW",
        "btag_bc_down_uncorrelated" : "BTAG_BC_UNCORR_DW",
        #"btag_l_UP"  : "BTAG_L_UP",
        "btag_l_up" : "BTAG_L_UP",
        "btag_l_up_correlated" : "BTAG_L_CORR_UP",
        "btag_l_up_uncorrelated" : "BTAG_L_UNCORR_UP",
        #"btag_l_DW"  : "BTAG_L_DW",
        "btag_l_down" : "BTAG_L_DW",
        "btag_l_down_correlated" : "BTAG_L_CORR_DW",
        "btag_l_down_uncorrelated" : "BTAG_L_UNCORR_DW",
        "PileupUp"  : "PILEUP_UP",
        "PileupDown": "PILEUP_DW",
        "PrefireUp" : "PREFIRE_UP",
        "PrefireDown" : "PREFIRE_DW",
        "ISRUp"  : "ISR_UP",
        "ISRDown": "ISR_DW",
        "FSRUp"  : "FSR_UP",
        "FSRDown": "FSR_DW",
        "FACTORUp"  : "FACTOR_UP",
        "FACTORDown": "FACTOR_DW",
        "RENORMUp"  : "RENORM_UP",
        "RENORMDown": "RENORM_DW",
        "RENORM_FACTOR_SAMEUp"  : "RENORM_FACTOR_UP",
        "RENORM_FACTOR_SAMEDown": "RENORM_FACTOR_DW",
        "RENORM_FACTOR_DIFFUp"  : "RENORM_UP_FACTOR_DW",
        "RENORM_FACTOR_DIFFDown": "RENORM_DW_FACTOR_UP",
        "mtopUP"  : "MTOP_UP",
        "mtopDOWN": "MTOP_DW",
        "hdampUP"  : "HDAMP_UP",
        "hdampDOWN": "HDAMP_DW",
        "ueUP"  : "UE_UP",
        "ueDOWN": "UE_DW",
        "JER_UP" : "JER_UP",
        "JER_DW" : "JER_DW",
        "JES_AbsoluteMPFBias_UP" : "JES_AbsoluteMPFBias_UP",
        "JES_AbsoluteMPFBias_DW" : "JES_AbsoluteMPFBias_DW",
        "JES_AbsoluteScale_UP" : "JES_AbsoluteScale_UP",
        "JES_AbsoluteScale_DW" : "JES_AbsoluteScale_DW",
        "JES_AbsoluteStat_UP" : "JES_AbsoluteStat_UP",
        "JES_AbsoluteStat_DW" : "JES_AbsoluteStat_DW",
        "JES_FlavorQCD_UP" : "JES_FlavorQCD_UP",
        "JES_FlavorQCD_DW" : "JES_FlavorQCD_DW",
        "JES_Fragmentation_UP" : "JES_Fragmentation_UP",
        "JES_Fragmentation_DW" : "JES_Fragmentation_DW",
        "JES_PileUpDataMC_UP" : "JES_PileUpDataMC_UP",
        "JES_PileUpDataMC_DW" : "JES_PileUpDataMC_DW",
        "JES_PileUpPtBB_UP" : "JES_PileUpPtBB_UP",
        "JES_PileUpPtBB_DW" : "JES_PileUpPtBB_DW",
        "JES_PileUpPtEC1_UP" : "JES_PileUpPtEC1_UP",
        "JES_PileUpPtEC1_DW" : "JES_PileUpPtEC1_DW",
        "JES_PileUpPtRef_UP" : "JES_PileUpPtRef_UP",
        "JES_PileUpPtRef_DW" : "JES_PileUpPtRef_DW",
        "JES_RelativeBal_UP" : "JES_RelativeBal_UP",
        "JES_RelativeBal_DW" : "JES_RelativeBal_DW",
        "JES_RelativeFSR_UP" : "JES_RelativeFSR_UP",
        "JES_RelativeFSR_DW" : "JES_RelativeFSR_DW",
        "JES_RelativeJEREC1_UP" : "JES_RelativeJEREC1_UP",
        "JES_RelativeJEREC1_DW" : "JES_RelativeJEREC1_DW",
        "JES_RelativePtBB_UP" : "JES_RelativePtBB_UP",
        "JES_RelativePtBB_DW" : "JES_RelativePtBB_DW",
        "JES_RelativePtEC1_UP" : "JES_RelativePtEC1_UP",
        "JES_RelativePtEC1_DW" : "JES_RelativePtEC1_DW",
        "JES_RelativeSample_UP" : "JES_RelativeSample_UP",
        "JES_RelativeSample_DW" : "JES_RelativeSample_DW",
        "JES_RelativeStatEC_UP" : "JES_RelativeStatEC_UP",
        "JES_RelativeStatEC_DW" : "JES_RelativeStatEC_DW",
        "JES_RelativeStatFSR_UP" : "JES_RelativeStatFSR_UP",
        "JES_RelativeStatFSR_DW" : "JES_RelativeStatFSR_DW",
        "JES_SinglePionECAL_UP" : "JES_SinglePionECAL_UP",
        "JES_SinglePionECAL_DW" : "JES_SinglePionECAL_DW",
        "JES_SinglePionHCAL_UP" : "JES_SinglePionHCAL_UP",
        "JES_SinglePionHCAL_DW" : "JES_SinglePionHCAL_DW",
        "JES_TimePtEta_UP" : "JES_TimePtEta_UP",
        "JES_TimePtEta_DW" : "JES_TimePtEta_DW",
        "JES_Total_UP" : "JES_Total_UP",
        "JES_Total_DW" : "JES_Total_DW",
        #"JES_UP" : "JES_UP",
        #"JES_DW" : "JES_DW",
        "MET_UP" : "MET_UP",
        "MET_DW" : "MET_DW",
        "mtop1665": "MTOP6GEV_DW",
        "mtop1695": "MTOP3GEV_DW",
        "mtop1715": "MTOP1GEV_DW",
        "mtop1735": "MTOP1GEV_UP",
        "mtop1755": "MTOP3GEV_UP",
        "mtop1785": "MTOP6GEV_UP",
            # A/H signal systematics
        "AH_RENORMUp"   : "AH_RENORM_UP",
        "AH_RENORMDown" : "AH_RENORM_DW",
        "AH_FACTORUp"   : "AH_FACTOR_UP",
        "AH_FACTORDown" : "AH_FACTOR_DW",
    }
    for year in ["2016APV", "2016", "2017", "2018"]
}
sys_to_name["2018"].pop("PrefireUp")
sys_to_name["2018"].pop("PrefireDown")



    ## dict for writing template names
template_sys_to_name = {
        ## format is name_in_analyzer : (extension_to_topology, only_applies_to_TT)
    year : {
        "nosys" : "",
        "JES_AbsoluteMPFBias_UP" : "JesAbsoluteMPFBiasUp",
        "JES_AbsoluteMPFBias_DW" : "JesAbsoluteMPFBiasDown",
        "JES_AbsoluteScale_UP" : "JesAbsoluteScaleUp",
        "JES_AbsoluteScale_DW" : "JesAbsoluteScaleDown",
        "JES_AbsoluteStat_UP" : "JesAbsoluteStatUp",
        "JES_AbsoluteStat_DW" : "JesAbsoluteStatDown",
        "JES_FlavorQCD_UP" : "JesFlavorQCDUp",
        "JES_FlavorQCD_DW" : "JesFlavorQCDDown",
        "JES_Fragmentation_UP" : "JesFragmentationUp",
        "JES_Fragmentation_DW" : "JesFragmentationDown",
        "JES_PileUpDataMC_UP" : "JesPileUpDataMCUp",
        "JES_PileUpDataMC_DW" : "JesPileUpDataMCDown",
        "JES_PileUpPtBB_UP" : "JesPileUpPtBBUp",
        "JES_PileUpPtBB_DW" : "JesPileUpPtBBDown",
        "JES_PileUpPtEC1_UP" : "JesPileUpPtEC1Up",
        "JES_PileUpPtEC1_DW" : "JesPileUpPtEC1Down",
        "JES_PileUpPtRef_UP" : "JesPileUpPtRefUp",
        "JES_PileUpPtRef_DW" : "JesPileUpPtRefDown",
        "JES_RelativeBal_UP" : "JesRelativeBalUp",
        "JES_RelativeBal_DW" : "JesRelativeBalDown",
        "JES_RelativeFSR_UP" : "JesRelativeFSRUp",
        "JES_RelativeFSR_DW" : "JesRelativeFSRDown",
        "JES_RelativeJEREC1_UP" : "JesRelativeJEREC1Up",
        "JES_RelativeJEREC1_DW" : "JesRelativeJEREC1Down",
        "JES_RelativePtBB_UP" : "JesRelativePtBBUp",
        "JES_RelativePtBB_DW" : "JesRelativePtBBDown",
        "JES_RelativePtEC1_UP" : "JesRelativePtEC1Up",
        "JES_RelativePtEC1_DW" : "JesRelativePtEC1Down",
        "JES_RelativeSample_UP" : "JesRelativeSampleUp",
        "JES_RelativeSample_DW" : "JesRelativeSampleDown",
        "JES_RelativeStatEC_UP" : "JesRelativeStatECUp",
        "JES_RelativeStatEC_DW" : "JesRelativeStatECDown",
        "JES_RelativeStatFSR_UP" : "JesRelativeStatFSRUp",
        "JES_RelativeStatFSR_DW" : "JesRelativeStatFSRDown",
        "JES_SinglePionECAL_UP" : "JesSinglePionECALUp",
        "JES_SinglePionECAL_DW" : "JesSinglePionECALDown",
        "JES_SinglePionHCAL_UP" : "JesSinglePionHCALUp",
        "JES_SinglePionHCAL_DW" : "JesSinglePionHCALDown",
        "JES_TimePtEta_UP" : "JesTimePtEtaUp",
        "JES_TimePtEta_DW" : "JesTimePtEtaDown",
        "JES_Total_UP" : "JesTotalUp",
        "JES_Total_DW" : "JesTotalDown",
        #"JES_UP" : "JesUp",
        #"JES_DW" : "JesDown",
        "JER_UP" : "JerUp",
        "JER_DW" : "JerDown",
        "MET_UP" : "UncMETUp",
        "MET_DW" : "UncMETDown",
        #"btag_bc_UP" : "btagUp",
        "btag_bc_up" : "btagUp",
        "btag_bc_up_correlated" : "btagCorrUp",
        "btag_bc_up_uncorrelated" : "btagUncorrUp",
        #"btag_bc_DW" : "btagDown",
        "btag_bc_down" : "btagDown",
        "btag_bc_down_correlated" : "btagCorrDown",
        "btag_bc_down_uncorrelated" : "btagUncorrDown",
        #"btag_l_UP" : "ltagUp",
        "btag_l_up" : "ltagUp",
        "btag_l_up_correlated" : "ltagCorrUp",
        "btag_l_up_uncorrelated" : "ltagUncorrUp",
        #"btag_l_DW" : "ltagDown",
        "btag_l_down" : "ltagDown",
        "btag_l_down_correlated" : "ltagCorrDown",
        "btag_l_down_uncorrelated" : "ltagUncorrDown",
        "Lep_RECOUp"   : "LEPsfrecoUp",
        "Lep_RECODown" : "LEPsfrecoDown",
        "Lep_TRIGUp"   : "LEPsftriggerUp",
        "Lep_TRIGDown" : "LEPsftriggerDown",
        "PileupUp"   : "pileupUp",
        "PileupDown" : "pileupDown",
        "PrefireUp"   : "prefireUp",
        "PrefireDown" : "prefireDown",
            # TT systematics
        "EWcorrUp" : "deltaQCDdeltaEW",
        #"EWcorrUp" : "EWcorr",
        "YukawaUp" : "ytUp",
        "YukawaDown" : "ytDown",
        "FSRUp"   : "psfsrUp",
        "FSRDown" : "psfsrDown",
        "ISRUp"   : "psisrUp",
        "ISRDown" : "psisrDown",
        "RENORMUp"   : "tt_uRUp",
        "RENORMDown" : "tt_uRDown",
        "FACTORUp"   : "tt_uFUp",
        "FACTORDown" : "tt_uFDown",
        #"RENORM_FACTOR_SAMEUp"   : "tt_uFuRUp",
        #"RENORM_FACTOR_SAMEDown" : "tt_uFuRDown",
        "ueUP"   : "uetuneUp",
        "ueDOWN" : "uetuneDown",
        "hdampUP"   : "hdampUp",
        "hdampDOWN" : "hdampDown",
        #"mtopUP"   : "TMassUp",
        #"mtopDOWN" : "TMassDown",
        "mtop1695" : "mtop3GeVDown",
        "mtop1755" : "mtop3GeVUp",
            # A/H signal systematics
        "AH_RENORMUp"   : "ah_uRUp",
        "AH_RENORMDown" : "ah_uRDown",
        "AH_FACTORUp"   : "ah_uFUp",
        "AH_FACTORDown" : "ah_uFDown",
    }
    for year in ["2016APV", "2016", "2017", "2018"]
}
template_sys_to_name["2018"].pop("PrefireUp")
template_sys_to_name["2018"].pop("PrefireDown")



# dict of systematic types with up/down variation names
sys_groups = {
    year : {
        #"nosys" : ["nosys"],
        "deltaQCDdeltaEW" : ["EWcorrUp", None],
        "YUKAWA" : ["YukawaUp", "YukawaDown"],
        "EWQCD_SHAPE" : ["shapeUp", "shapeDown"],
        "LEP_RECO" : ["Lep_RECOUp", "Lep_RECODown"],
        "LEP_TRIG" : ["Lep_TRIGUp", "Lep_TRIGDown"],
        #"BTAG_BC" : ["btag_bc_up", "btag_bc_down"],
        "BTAG_BC_CORR" : ["btag_bc_up_correlated", "btag_bc_down_correlated"],
        "BTAG_BC_UNCORR" : ["btag_bc_up_uncorrelated", "btag_bc_down_uncorrelated"],
        #"BTAG_L" : ["btag_l_up", "btag_l_down"],
        "BTAG_L_CORR" : ["btag_l_up_correlated", "btag_l_down_correlated"],
        "BTAG_L_UNCORR" : ["btag_l_up_uncorrelated", "btag_l_down_uncorrelated"],
        "PILEUP" : ["PileupUp", "PileupDown"],
        "PREFIRE" : ["PrefireUp", "PrefireDown"],
        "ISR" : ["ISRUp", "ISRDown"],
        "FSR" : ["FSRUp", "FSRDown"],
        "FACTOR" : ["FACTORUp", "FACTORDown"],
        "RENORM" : ["RENORMUp", "RENORMDown"],
        #"RENORM_FACTOR" : ["RENORM_FACTOR_SAMEUp", "RENORM_FACTOR_SAMEDown"],
        #"RENORM_UP_FACTOR_DW" : ["RENORM_FACTOR_DIFFUp", None],
        #"RENORM_DW_FACTOR_UP" : ["RENORM_FACTOR_DIFFDown", None],        
        #"MTOP" : ["mtopUP", "mtopDOWN"],
        "HDAMP" : ["hdampUP", "hdampDOWN"],
        "UE" : ["ueUP", "ueDOWN"],
        "JER" : ["JER_UP", "JER_DW"],
        "JES_AbsoluteMPFBias" : ["JES_AbsoluteMPFBias_UP", "JES_AbsoluteMPFBias_DW"],
        "JES_AbsoluteScale" : ["JES_AbsoluteScale_UP", "JES_AbsoluteScale_DW"],
        "JES_AbsoluteStat" : ["JES_AbsoluteStat_UP", "JES_AbsoluteStat_DW"],
        "JES_FlavorQCD" : ["JES_FlavorQCD_UP", "JES_FlavorQCD_DW"],
        "JES_Fragmentation" : ["JES_Fragmentation_UP", "JES_Fragmentation_DW"],
        "JES_PileUpDataMC" : ["JES_PileUpDataMC_UP", "JES_PileUpDataMC_DW"],
        "JES_PileUpPtBB" : ["JES_PileUpPtBB_UP", "JES_PileUpPtBB_DW"],
        "JES_PileUpPtEC1" : ["JES_PileUpPtEC1_UP", "JES_PileUpPtEC1_DW"],
        "JES_PileUpPtRef" : ["JES_PileUpPtRef_UP", "JES_PileUpPtRef_DW"],
        "JES_RelativeBal" : ["JES_RelativeBal_UP", "JES_RelativeBal_DW"],
        "JES_RelativeFSR" : ["JES_RelativeFSR_UP", "JES_RelativeFSR_DW"],
        "JES_RelativeJEREC1" : ["JES_RelativeJEREC1_UP", "JES_RelativeJEREC1_DW"],
        "JES_RelativePtBB" : ["JES_RelativePtBB_UP", "JES_RelativePtBB_DW"],
        "JES_RelativePtEC1" : ["JES_RelativePtEC1_UP", "JES_RelativePtEC1_DW"],
        "JES_RelativeSample" : ["JES_RelativeSample_UP", "JES_RelativeSample_DW"],
        "JES_RelativeStatEC" : ["JES_RelativeStatEC_UP", "JES_RelativeStatEC_DW"],
        "JES_RelativeStatFSR" : ["JES_RelativeStatFSR_UP", "JES_RelativeStatFSR_DW"],
        "JES_SinglePionECAL" : ["JES_SinglePionECAL_UP", "JES_SinglePionECAL_DW"],
        "JES_SinglePionHCAL" : ["JES_SinglePionHCAL_UP", "JES_SinglePionHCAL_DW"],
        "JES_TimePtEta" : ["JES_TimePtEta_UP", "JES_TimePtEta_DW"],
        "JES_Total" : ["JES_Total_UP", "JES_Total_DW"],
        "MET" : ["MET_UP", "MET_DW"],
        "MTOP1GEV" : ["mtop1735", "mtop1715"],
        "MTOP3GEV" : ["mtop1755", "mtop1695"],
        "MTOP6GEV" : ["mtop1785", "mtop1665"],
            # A/H signal systematics
        "AH_RENORM" : ["AH_RENORMUp", "AH_RENORMDown"],
        "AH_FACTOR" : ["AH_FACTORUp", "AH_FACTORDown"],
    }
    for year in ["2016APV", "2016", "2017", "2018"]
}
sys_groups["2018"].pop("PREFIRE")


