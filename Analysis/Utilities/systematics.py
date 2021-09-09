"""
List of systematics
"""

#    ## individual types
#lep_sys = [
#    "Lep_RECOUp", "Lep_RECODown",
#    "Lep_TRIGUp", "Lep_TRIGDown"
#]
#
#btag_sys = [
#    "btag_bc_UP", "btag_bc_DW",
#    "btag_l_UP", "btag_l_DW",
#]
#
#jet_sys = [
#    "JER_UP", "JER_DW",
#    "JES_UP", "JES_DW",
#    "MET_UP", "MET_DW",
#]

ind_jec_sys = {
    year : {
        #"JES_AbsoluteFlavMap",
        "JES_AbsoluteMPFBias",
        "JES_AbsoluteScale",
        "JES_AbsoluteStat",
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
        "JES_FlavorQCD",
        #"JES_FlavorZJet",
        "JES_Fragmentation",
        "JES_PileUpDataMC",
        #"JES_PileUpEnvelope",
        #"JES_PileUpMuZero",
        "JES_PileUpPtBB",
        "JES_PileUpPtEC1",
        #"JES_PileUpPtEC2",
        #"JES_PileUpPtHF",
        "JES_PileUpPtRef",
        "JES_RelativeBal",
        "JES_RelativeFSR",
        "JES_RelativeJEREC1",
        #"JES_RelativeJEREC2",
        #"JES_RelativeJERHF",
        "JES_RelativePtBB",
        "JES_RelativePtEC1",
        #"JES_RelativePtEC2",
        #"JES_RelativePtHF",
        "JES_RelativeSample",
        "JES_RelativeStatEC",
        "JES_RelativeStatFSR",
        #"JES_RelativeStatHF",
        "JES_SinglePionECAL",
        "JES_SinglePionHCAL",
        #"JES_SubTotalAbsolute",
        #"JES_SubTotalMC",
        #"JES_SubTotalPileUp",
        #"JES_SubTotalPt",
        #"JES_SubTotalRelative",
        #"JES_SubTotalScale",
        "JES_TimePtEta",
        "JES_Total",
        #"JES_TotalNoFlavor",
        #"JES_TotalNoFlavorNoTime",
        #"JES_TotalNoTime",
    }
    for year in ["2016APV", "2016", "2017", "2018"]
}
#ind_jec_sys["2016APV"].update({
#    "JES_TimeRunBCD",
#    "JES_TimeRunEF",
#    "JES_TimeRunGH",
#})
#ind_jec_sys["2016"].update({
#    "JES_TimeRunBCD",
#    "JES_TimeRunEF",
#    "JES_TimeRunGH",
#})
#ind_jec_sys["2017"].update({
#    "JES_TimeRunB",
#    "JES_TimeRunC",
#    "JES_TimeRunDE",
#    "JES_TimeRunF",
#})
#ind_jec_sys["2018"].update({
#    "JES_AbsoluteSample",
#    "JES_TimeRunA",
#    "JES_TimeRunB",
#    "JES_TimeRunC",
#    "JES_TimeRunD",
#})



#misc_sys = [
#    "PileupUp", "PileupDown",
#]
#
#prefire_sys = ["PrefireUp", "PrefireDown"]
#
#ps_sys = [
#    "ISRUp", "ISRDown",
#    "FSRUp", "FSRDown",
#]
#
#lhe_sys = [
#    "FACTORUp", "FACTORDown",
#    "RENORMUp", "RENORMDown",
#    "RENORM_FACTOR_SAMEUp", "RENORM_FACTOR_SAMEDown",
#    "RENORM_FACTOR_DIFFUp", "RENORM_FACTOR_DIFFDown",
#]
#    ## groups
#reweight_sys = {
#    "2016APV" : lep_sys + btag_sys + misc_sys + prefire_sys,
#    "2016" : lep_sys + btag_sys + misc_sys + prefire_sys,
#    "2017" : lep_sys + btag_sys + misc_sys + prefire_sys,
#    "2018" : lep_sys + btag_sys + misc_sys
#}

#event_sys = {
#    "2016APV" : jet_sys,
#    "2016" : jet_sys,
#    "2017" : jet_sys,
#    "2018" : jet_sys
#}

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
        "BTAG_BC_UP" : "btag_bc_UP",
        "BTAG_BC_DW" : "btag_bc_DW",
        "BTAG_L_UP" : "btag_l_UP",
        "BTAG_L_DW" : "btag_l_DW",
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


# dict of name_I_want : name_in_code
    # systematics available to be applied to signal (reweighting)
signal_reweight_opts = {
    year : {
        "LEP_RECO_UP" : "Lep_RECOUp",
        "LEP_RECO_DW" : "Lep_RECODown",
        "LEP_TRIG_UP" : "Lep_TRIGUp",
        "LEP_TRIG_DW" : "Lep_TRIGDown",
        "BTAG_BC_UP" : "btag_bc_UP",
        "BTAG_BC_DW" : "btag_bc_DW",
        "BTAG_L_UP" : "btag_l_UP",
        "BTAG_L_DW" : "btag_l_DW",
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
    "2016APV" : {
        "nosys" : "nosys",
        "Lep_RECOUp"   : "LEP_RECO_UP",
        "Lep_RECODown" : "LEP_RECO_DW",
        "Lep_TRIGUp"   : "LEP_TRIG_UP",
        "Lep_TRIGDown" : "LEP_TRIG_DW",
        "btag_bc_UP" : "BTAG_BC_UP",
        "btag_bc_DW" : "BTAG_BC_DW",
        "btag_l_UP"  : "BTAG_L_UP",
        "btag_l_DW"  : "BTAG_L_DW",
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
        "JES_UP" : "JES_UP",
        "JES_DW" : "JES_DW",
        "MET_UP" : "MET_UP",
        "MET_DW" : "MET_DW",
        "fsrUP"  : "sampleFSR_UP",
        "fsrDOWN": "sampleFSR_DW",
        "isrUP"  : "sampleISR_UP",
        "isrDOWN": "sampleISR_DW",
        "mtop1665": "MTOP6GeV_DW",
        "mtop1695": "MTOP3GeV_DW",
        "mtop1715": "MTOP1GeV_DW",
        "mtop1735": "MTOP1GeV_UP",
        "mtop1755": "MTOP3GeV_UP",
        "mtop1785": "MTOP6GeV_UP",
    },
    "2016" : {
        "nosys" : "nosys",
        "Lep_RECOUp"   : "LEP_RECO_UP",
        "Lep_RECODown" : "LEP_RECO_DW",
        "Lep_TRIGUp"   : "LEP_TRIG_UP",
        "Lep_TRIGDown" : "LEP_TRIG_DW",
        "btag_bc_UP" : "BTAG_BC_UP",
        "btag_bc_DW" : "BTAG_BC_DW",
        "btag_l_UP"  : "BTAG_L_UP",
        "btag_l_DW"  : "BTAG_L_DW",
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
        "JES_UP" : "JES_UP",
        "JES_DW" : "JES_DW",
        "MET_UP" : "MET_UP",
        "MET_DW" : "MET_DW",
        "fsrUP"  : "sampleFSR_UP",
        "fsrDOWN": "sampleFSR_DW",
        "isrUP"  : "sampleISR_UP",
        "isrDOWN": "sampleISR_DW",
        "mtop1665": "MTOP6GeV_DW",
        "mtop1695": "MTOP3GeV_DW",
        "mtop1715": "MTOP1GeV_DW",
        "mtop1735": "MTOP1GeV_UP",
        "mtop1755": "MTOP3GeV_UP",
        "mtop1785": "MTOP6GeV_UP",
    },
    "2017" : {
        "nosys" : "nosys",
        "Lep_RECOUp"   : "LEP_RECO_UP",
        "Lep_RECODown" : "LEP_RECO_DW",
        "Lep_TRIGUp"   : "LEP_TRIG_UP",
        "Lep_TRIGDown" : "LEP_TRIG_DW",
        "btag_bc_UP" : "BTAG_BC_UP",
        "btag_bc_DW" : "BTAG_BC_DW",
        "btag_l_UP"  : "BTAG_L_UP",
        "btag_l_DW"  : "BTAG_L_DW",
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
        "JES_UP" : "JES_UP",
        "JES_DW" : "JES_DW",
        "MET_UP" : "MET_UP",
        "MET_DW" : "MET_DW",
        "mtop1665": "MTOP6GeV_DW",
        "mtop1695": "MTOP3GeV_DW",
        "mtop1715": "MTOP1GeV_DW",
        "mtop1735": "MTOP1GeV_UP",
        "mtop1755": "MTOP3GeV_UP",
        "mtop1785": "MTOP6GeV_UP",
    },
    "2018" : {
        "nosys" : "nosys",
        "Lep_RECOUp"   : "LEP_RECO_UP",
        "Lep_RECODown" : "LEP_RECO_DW",
        "Lep_TRIGUp"   : "LEP_TRIG_UP",
        "Lep_TRIGDown" : "LEP_TRIG_DW",
        "btag_bc_UP" : "BTAG_BC_UP",
        "btag_bc_DW" : "BTAG_BC_DW",
        "btag_l_UP"  : "BTAG_L_UP",
        "btag_l_DW"  : "BTAG_L_DW",
        "PileupUp"  : "PILEUP_UP",
        "PileupDown": "PILEUP_DW",
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
        "JES_UP" : "JES_UP",
        "JES_DW" : "JES_DW",
        "MET_UP" : "MET_UP",
        "MET_DW" : "MET_DW",
        "mtop1665": "MTOP6GeV_DW",
        "mtop1695": "MTOP3GeV_DW",
        "mtop1715": "MTOP1GeV_DW",
        "mtop1735": "MTOP1GeV_UP",
        "mtop1755": "MTOP3GeV_UP",
        "mtop1785": "MTOP6GeV_UP",
    },
}



    ## dict for writing template names
template_sys_to_name = {
        ## format is name_in_analyzer : (extension_to_topology, only_applies_to_TT)
    #year : {
    #    "nosys" : ("", False),
    #    "JES_UP" : ("CMS_scale_j_13TeVUp", False),
    #    "JES_DW" : ("CMS_scale_j_13TeVDown", False),
    #    "JER_UP" : ("CMS_res_j_13TeVUp", False),
    #    "JER_DW" : ("CMS_res_j_13TeVDown", False),
    #    "MET_UP" : ("CMS_METunclustered_13TeVUp", False),
    #    "MET_DW" : ("CMS_METunclustered_13TeVDown", False),
    #    "btag_bc_UP" : ("CMS_eff_b_13TeVUp", False),
    #    "btag_bc_DW" : ("CMS_eff_b_13TeVDown", False),
    #    "btag_l_UP" : ("CMS_fake_b_13TeVUp", False),
    #    "btag_l_DW" : ("CMS_fake_b_13TeVDown", False),
    #    "Lep_RECOUp" : ("CMS_eff_reco_LEPUp", False),
    #    "Lep_RECODown" : ("CMS_eff_reco_LEPDown", False),
    #    "Lep_TRIGUp" : ("CMS_eff_trigger_LEPUp", False),
    #    "Lep_TRIGDown" : ("CMS_eff_trigger_LEPDown", False),
    #    "PileupUp" : ("CMS_pileupUp", False),
    #    "PileupDown" : ("CMS_pileupDown", False),
    #    "PrefireUp" : ("CMS_prefireUp", False),
    #    "PrefireDown" : ("CMS_prefireDown", False), 
    #        # TT systematics
    #    "FSRUp" : ("FSR_TTUp", True),
    #    "FSRDown" : ("FSR_TTDown", True),
    #    "ISRUp" : ("ISR_TTUp", True),
    #    "ISRDown" : ("ISR_TTDown", True),
    #    "RENORMUp" : ("QCDscaleMERenorm_TTUp", True),
    #    "RENORMDown" : ("QCDscaleMERenorm_TTDown", True),
    #    "FACTORUp" : ("QCDscaleMEFactor_TTUp", True),
    #    "FACTORDown" : ("QCDscaleMEFactor_TTDown", True),
    #    "RENORM_FACTOR_SAMEUp" : ("QCDscaleMERenormFactor_TTUp", True),
    #    "RENORM_FACTOR_SAMEDown" : ("QCDscaleMERenormFactor_TTDown", True),
    #    "ueUP" : ("UE_TTUp", True),
    #    "ueDOWN" : ("UE_TTDown", True),
    #    "hdampUP" : ("Hdamp_TTUp", True),
    #    "hdampDOWN" : ("Hdamp_TTDown", True),
    #    #"mtopUP" : ("TMassUp", True),
    #    #"mtopDOWN" : ("TMassDown", True),
    #    "mtop1695" : ("TMassDown", True),
    #    "mtop1755" : ("TMassUp", True),
    #    #"mtop1695" : ("TMass3GeVDown", True),
    #    #"mtop1755" : ("TMass3GeVUp", True),
    #    #"mtop1695" : ("TMassDown", True),
    #    #"mtop1755" : ("TMassUp", True),
    #}
    year : {
        "nosys" : ("", False),
        "JES_UP" : ("JesUp", False),
        "JES_DW" : ("JesDown", False),
        "JER_UP" : ("JerUp", False),
        "JER_DW" : ("JerDown", False),
        "MET_UP" : ("UncMETUp", False),
        "MET_DW" : ("UncMETDown", False),
        "btag_bc_UP" : ("btagUp", False),
        "btag_bc_DW" : ("btagDown", False),
        "btag_l_UP" : ("ltagUp", False),
        "btag_l_DW" : ("ltagDown", False),
        "Lep_RECOUp" : ("LEPsfrecoUp", False),
        "Lep_RECODown" : ("LEPsfrecoDown", False),
        "Lep_TRIGUp" : ("LEPsftriggerUp", False),
        "Lep_TRIGDown" : ("LEPsftriggerDown", False),
        "PileupUp" : ("pileupUp", False),
        "PileupDown" : ("pileupDown", False),
        "PrefireUp" : ("prefireUp", False),
        "PrefireDown" : ("prefireDown", False), 
            # TT systematics
        "FSRUp" : ("psfsrUp", True),
        "FSRDown" : ("psfsrDown", True),
        "ISRUp" : ("psisrUp", True),
        "ISRDown" : ("psisrDown", True),
        "RENORMUp" : ("tt_uRUp", True),
        "RENORMDown" : ("tt_uRDown", True),
        "FACTORUp" : ("tt_uFUp", True),
        "FACTORDown" : ("tt_uFDown", True),
        #"RENORM_FACTOR_SAMEUp" : ("tt_uFuRUp", True),
        #"RENORM_FACTOR_SAMEDown" : ("tt_uFuRDown", True),
        "ueUP" : ("uetuneUp", True),
        "ueDOWN" : ("uetuneDown", True),
        "hdampUP" : ("hdampUp", True),
        "hdampDOWN" : ("hdampDown", True),
        #"mtopUP" : ("TMassUp", True),
        #"mtopDOWN" : ("TMassDown", True),
        "mtop1695" : ("mtop3GeVDown", True),
        "mtop1755" : ("mtop3GeVUp", True),
            # A/H signal systematics
        "AH_RENORMUp" : ("ah_uRUp", True),
        "AH_RENORMDown" : ("ah_uRDown", True),
        "AH_FACTORUp" : ("ah_uFUp", True),
        "AH_FACTORDown" : ("ah_uFDown", True),
    }
    for year in ["2016APV", "2016", "2017", "2018"]
}
template_sys_to_name["2018"].pop("PrefireUp")
template_sys_to_name["2018"].pop("PrefireDown")

