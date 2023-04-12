from coffea.lookup_tools.dense_lookup import dense_lookup
from pdb import set_trace
import numpy as np

class LeptonSF(object):
    def __init__(self, input_SFs = None, debug = False):

        self.input_SFs = input_SFs
        self.debug = debug

        self.schema_ = {
            "Muons": {
                "central" : ["ID_Central", "ISO_Central", "TRIG_Central", "RECO_Central"],
                    # variations of ID
                "IDtotUp" : ["ID_Error_totUp", "ISO_Central", "TRIG_Central", "RECO_Central"],
                "IDtotDown" : ["ID_Error_totDown", "ISO_Central", "TRIG_Central", "RECO_Central"],
                "IDstatUp" : ["ID_Error_statUp", "ISO_Central", "TRIG_Central", "RECO_Central"],
                "IDstatDown" : ["ID_Error_statDown", "ISO_Central", "TRIG_Central", "RECO_Central"],
                "IDsystUp" : ["ID_Error_systUp", "ISO_Central", "TRIG_Central", "RECO_Central"],
                "IDsystDown" : ["ID_Error_systDown", "ISO_Central", "TRIG_Central", "RECO_Central"],
                    # variations of ISO
                "ISOtotUp" : ["ID_Central", "ISO_Error_totUp", "TRIG_Central", "RECO_Central"],
                "ISOtotDown" : ["ID_Central", "ISO_Error_totDown", "TRIG_Central", "RECO_Central"],
                "ISOstatUp" : ["ID_Central", "ISO_Error_statUp", "TRIG_Central", "RECO_Central"],
                "ISOstatDown" : ["ID_Central", "ISO_Error_statDown", "TRIG_Central", "RECO_Central"],
                "ISOsystUp" : ["ID_Central", "ISO_Error_systUp", "TRIG_Central", "RECO_Central"],
                "ISOsystDown" : ["ID_Central", "ISO_Error_systDown", "TRIG_Central", "RECO_Central"],
                    # variations of TRIG
                "TRIGtotUp" : ["ID_Central", "ISO_Central", "TRIG_Error_totUp", "RECO_Central"],
                "TRIGtotDown" : ["ID_Central", "ISO_Central", "TRIG_Error_totDown", "RECO_Central"],
                "TRIGstatUp" : ["ID_Central", "ISO_Central", "TRIG_Error_statUp", "RECO_Central"],
                "TRIGstatDown" : ["ID_Central", "ISO_Central", "TRIG_Error_statDown", "RECO_Central"],
                "TRIGsystUp" : ["ID_Central", "ISO_Central", "TRIG_Error_systUp", "RECO_Central"],
                "TRIGsystDown" : ["ID_Central", "ISO_Central", "TRIG_Error_systDown", "RECO_Central"],
                    # variations of RECO
                "RECOtotUp" : ["ID_Central", "ISO_Central", "TRIG_Central", "RECO_Error_totUp"],
                "RECOtotDown" : ["ID_Central", "ISO_Central", "TRIG_Central", "RECO_Error_totDown"],
            },
            "Electrons" : {
                "central" : ["ID_Central", "TRIG_Central", "RECO_Central"],
                    # variations of ID
                "IDtotUp" : ["ID_Error_totUp", "TRIG_Central", "RECO_Central"],
                "IDtotDown" : ["ID_Error_totDown", "TRIG_Central", "RECO_Central"],
                    # variations of TRIG
                "TRIGtotUp" : ["ID_Central", "TRIG_Error_totUp", "RECO_Central"],
                "TRIGtotDown" : ["ID_Central", "TRIG_Error_totDown", "RECO_Central"],
                    # variations of RECO
                "RECOtotUp" : ["ID_Central", "TRIG_Central", "RECO_Error_totUp"],
                "RECOtotDown" : ["ID_Central", "TRIG_Central", "RECO_Error_totDown"],
            }
        }

        self.eta_edges = {"Muons" : np.around(np.linspace(-2.5, 2.5, 101), decimals=3), "Electrons" : np.around(np.linspace(-2.5, 2.5, 101), decimals=3)}
        self.pt_edges = {"Muons" : np.around(np.linspace(0., 130., 27), decimals=1), "Electrons" : np.around(np.linspace(0., 200., 41), decimals=1)}
        self.eta_centers = {lep : np.array([(self.eta_edges[lep][idx]+self.eta_edges[lep][idx+1])/2 for idx in range(len(self.eta_edges[lep])-1)]) for lep in self.eta_edges.keys()}
        self.pt_centers = {lep : np.array([(self.pt_edges[lep][idx]+self.pt_edges[lep][idx+1])/2 for idx in range(len(self.pt_edges[lep])-1)]) for lep in self.pt_edges.keys()}

        self.indivSF_sources_ = self.indiv_variations_()
        self.combSF_sources_ = self.combined_variations_()


    def indiv_variations_(self):
        #set_trace()
        indiv_SF_sources = {"Muons" : {}, "Electrons" : {}}
        for lepflav in indiv_SF_sources.keys():
            pt_vals, eta_vals = np.meshgrid(self.pt_centers[lepflav], self.eta_centers[lepflav])
            #if lepflav == "Electrons": set_trace()
            #if lepflav == "Muons": set_trace()
            for sf_type in self.input_SFs[lepflav].keys():
                isAbsEta = self.input_SFs[lepflav][sf_type]["isAbsEta"]
                if "eta_ranges" in self.input_SFs[lepflav][sf_type].keys():
                    eta_ranges = self.input_SFs[lepflav][sf_type]["eta_ranges"]
                    for sf_var in self.input_SFs[lepflav][sf_type].keys():
                        if ((sf_var == "eta_ranges") or (sf_var == "isAbsEta")): continue
                        if self.debug: print(lepflav, sf_type, sf_var)
                        if sf_var == "Central":
                            cen_wts = np.ones(pt_vals.shape)
                        elif "Error" in sf_var:
                            errup_wts, errdw_wts = np.ones(pt_vals.shape), np.ones(pt_vals.shape)
                        else:
                            set_trace()
    
                        for idx, eta_range in enumerate(eta_ranges):
                            mask = (np.abs(eta_vals) >= eta_range[0]) & (np.abs(eta_vals) < eta_range[1]) if isAbsEta else (eta_vals >= eta_range[0]) & (eta_vals < eta_range[1]) # find inds that are within given eta range
                            if not np.any(mask): continue # no values fall within eta range
                            if sf_var == "Central":
                                if self.input_SFs[lepflav][sf_type][sf_var][f"eta_bin{idx}"]._values.size > 0: # check that the dense_lookup has values
                                    cen_wts[mask] = self.input_SFs[lepflav][sf_type][sf_var][f"eta_bin{idx}"](pt_vals[mask])
                            else:
                                if (self.input_SFs[lepflav][sf_type]["Central"][f"eta_bin{idx}"]._values.size > 0) & (self.input_SFs[lepflav][sf_type][sf_var][f"eta_bin{idx}"]._values.size > 0): # check that the dense_lookup has values
                                    errup_wts[mask] = self.input_SFs[lepflav][sf_type]["Central"][f"eta_bin{idx}"](pt_vals[mask]) + self.input_SFs[lepflav][sf_type][sf_var][f"eta_bin{idx}"](pt_vals[mask])
                                    errdw_wts[mask] = self.input_SFs[lepflav][sf_type]["Central"][f"eta_bin{idx}"](pt_vals[mask]) - self.input_SFs[lepflav][sf_type][sf_var][f"eta_bin{idx}"](pt_vals[mask])
                        if sf_var == "Central":
                            indiv_SF_sources[lepflav][f"{sf_type}_{sf_var}"] = dense_lookup(cen_wts, (self.eta_edges[lepflav], self.pt_edges[lepflav]))
                        else:
                            indiv_SF_sources[lepflav][f"{sf_type}_{sf_var}Up"]   = dense_lookup(errup_wts, (self.eta_edges[lepflav], self.pt_edges[lepflav]))
                            indiv_SF_sources[lepflav][f"{sf_type}_{sf_var}Down"] = dense_lookup(errdw_wts, (self.eta_edges[lepflav], self.pt_edges[lepflav]))
                else:
                    for sf_var in self.input_SFs[lepflav][sf_type].keys():
                        if sf_var == "isAbsEta": continue
                        if self.debug: print(lepflav, sf_type, sf_var)
                        if sf_var == "Central":
                            if self.input_SFs[lepflav][sf_type][sf_var]._dimension == 1:
                                indiv_SF_sources[lepflav][f"{sf_type}_{sf_var}"] = dense_lookup(self.input_SFs[lepflav][sf_type][sf_var](np.abs(eta_vals)), (self.eta_edges[lepflav], self.pt_edges[lepflav])) if isAbsEta \
                                    else dense_lookup(self.input_SFs[lepflav][sf_type][sf_var](eta_vals), (self.eta_edges[lepflav], self.pt_edges[lepflav]))
                            elif self.input_SFs[lepflav][sf_type][sf_var]._dimension == 2:
                                indiv_SF_sources[lepflav][f"{sf_type}_{sf_var}"] = dense_lookup(self.input_SFs[lepflav][sf_type][sf_var](np.abs(eta_vals), pt_vals), (self.eta_edges[lepflav], self.pt_edges[lepflav])) if isAbsEta \
                                    else dense_lookup(self.input_SFs[lepflav][sf_type][sf_var](eta_vals, pt_vals), (self.eta_edges[lepflav], self.pt_edges[lepflav]))
                            else:
                                raise ValueError("Only 1D or 2D scale factors are supported!")
                        else:
                            if self.input_SFs[lepflav][sf_type][sf_var]._dimension == 1:
                                indiv_SF_sources[lepflav][f"{sf_type}_{sf_var}Up"]   = dense_lookup(self.input_SFs[lepflav][sf_type]["Central"](np.abs(eta_vals)) + self.input_SFs[lepflav][sf_type][sf_var](np.abs(eta_vals)), (self.eta_edges[lepflav], self.pt_edges[lepflav])) if isAbsEta \
                                    else dense_lookup(self.input_SFs[lepflav][sf_type]["Central"](eta_vals) + self.input_SFs[lepflav][sf_type][sf_var](eta_vals), (self.eta_edges[lepflav], self.pt_edges[lepflav]))
                                indiv_SF_sources[lepflav][f"{sf_type}_{sf_var}Down"] = dense_lookup(self.input_SFs[lepflav][sf_type]["Central"](np.abs(eta_vals)) - self.input_SFs[lepflav][sf_type][sf_var](np.abs(eta_vals)), (self.eta_edges[lepflav], self.pt_edges[lepflav])) if isAbsEta \
                                    else dense_lookup(self.input_SFs[lepflav][sf_type]["Central"](eta_vals) - self.input_SFs[lepflav][sf_type][sf_var](eta_vals), (self.eta_edges[lepflav], self.pt_edges[lepflav]))
                            elif self.input_SFs[lepflav][sf_type][sf_var]._dimension == 2:
                                indiv_SF_sources[lepflav][f"{sf_type}_{sf_var}Up"]   = dense_lookup(self.input_SFs[lepflav][sf_type]["Central"](np.abs(eta_vals), pt_vals) + self.input_SFs[lepflav][sf_type][sf_var](np.abs(eta_vals), pt_vals), (self.eta_edges[lepflav], self.pt_edges[lepflav])) if isAbsEta \
                                    else dense_lookup(self.input_SFs[lepflav][sf_type]["Central"](eta_vals, pt_vals) + self.input_SFs[lepflav][sf_type][sf_var](eta_vals, pt_vals), (self.eta_edges[lepflav], self.pt_edges[lepflav]))
                                indiv_SF_sources[lepflav][f"{sf_type}_{sf_var}Down"] = dense_lookup(self.input_SFs[lepflav][sf_type]["Central"](np.abs(eta_vals), pt_vals) - self.input_SFs[lepflav][sf_type][sf_var](np.abs(eta_vals), pt_vals), (self.eta_edges[lepflav], self.pt_edges[lepflav])) if isAbsEta \
                                    else dense_lookup(self.input_SFs[lepflav][sf_type]["Central"](eta_vals, pt_vals) - self.input_SFs[lepflav][sf_type][sf_var](eta_vals, pt_vals), (self.eta_edges[lepflav], self.pt_edges[lepflav]))
                            else:
                                raise ValueError("Only 1D or 2D scale factors are supported!")
    
        #set_trace()
        return indiv_SF_sources


    def combined_variations_(self):
        #set_trace()
        """
        Compute SF values corresponding to the key, value pairs in self.schema_ and save them as a dense_lookup object.
        """
        combSF_sources = {"Muons" : {}, "Electrons" : {}}
        for lepflav in combSF_sources.keys():
            for key, sf_list in self.schema_[lepflav].items():
                combSF_sources[lepflav][key] = dense_lookup(np.prod(np.dstack([self.indivSF_sources_[lepflav][sf_source]._values for sf_source in sf_list]), axis=2), (self.eta_edges[lepflav], self.pt_edges[lepflav]))
        return combSF_sources


    def get_scale_factor(self, pt, eta, tight_lep_mask, leptype, sysnames):
        #set_trace()
        output_SFs = {}
        for sys, dl in self.combSF_sources_[leptype].items():
            if sys not in sysnames+["central"]: continue
            evt_wts = np.ones(len(tight_lep_mask))
            evt_wts[tight_lep_mask] = dl(eta, pt)
            output_SFs[sys] = np.copy(evt_wts)

        return output_SFs
