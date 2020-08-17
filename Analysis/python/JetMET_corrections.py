'''
This script is supposed to be the same as https://github.com/CoffeaTeam/coffea/blob/master/coffea/jetmet_tools/JetTransformer.py
but handle up and down variations to jet and met separately for easier systematics computations
'''

from coffea.jetmet_tools.FactorizedJetCorrector import FactorizedJetCorrector
from coffea.jetmet_tools.JetResolution import JetResolution
from coffea.jetmet_tools.JetResolutionScaleFactor import JetResolutionScaleFactor
from coffea.jetmet_tools.JetCorrectionUncertainty import JetCorrectionUncertainty

from coffea.analysis_objects.JaggedCandidateArray import JaggedCandidateArray

import numpy as np
from uproot_methods.classes.TLorentzVector import (
    TLorentzVectorArray,
    ArrayMethods,
    PtEtaPhiMassArrayMethods
)
from copy import deepcopy
from pdb import set_trace
import warnings

_signature_map = {'JetPt': 'pt',
                  'JetEta': 'eta',
                  'Rho': 'rho',
                  'JetA': 'area'
                 }


def _update_jet_ptm(corr, jet, fromRaw=False):
    """
    This is a hack to update the jet pt and jet mass in place
    as we apply corrections and smearings.
    """
    if fromRaw:
        jet._content._contents['__fast_pt'] = corr * jet.ptRaw.content
        jet._content._contents['__fast_mass'] = corr * jet.massRaw.content
    else:
        jet._content._contents['__fast_pt'] = corr * jet.pt.content
        jet._content._contents['__fast_mass'] = corr * jet.mass.content

def _update_met_ptphi(met, pt, phi):
    """
    This is a hack to update the jet pt and jet mass in place
    as we apply corrections and smearings.
    """
    met._content._contents['__fast_pt'] = pt
    met._content._contents['__fast_phi'] = phi




def transform(corrections, jet, met, forceStochastic=False, shift=None):
    """
    This is the main entry point for JetTransformer and acts on arrays of jet data in-place.
    **precondition** : jet is a JaggedCandidateArray with additional attributes
        - 'ptRaw'
        - 'massRaw'
        - 'ptGenJet' if using hybrid JER
    You can call transform like this::
        xformer.transform(jets)
    **postcondition** : jet.pt, jet.mass, jet.p4 are updated to represent the corrected jet based on the input correction set.
    """
    if not isinstance(jet, JaggedCandidateArray):
        raise Exception('Input data must be a JaggedCandidateArray!')
    if ('ptRaw' not in jet.columns or 'massRaw' not in jet.columns):
        raise Exception('Input JaggedCandidateArray must have "ptRaw" & "massRaw"!')
    if ('ptGenJet' not in jet.columns):
        warnings.warn('Input JaggedCandidateArray must have "ptGenJet" in order to apply hybrid JER smearing method. Stochastic smearing will be applied.')
        forceStochastic = True

    if met is not None:
        if 'p4' not in met.columns:
            raise Exception('Input met must have a p4 column!')
        if not (isinstance(met['p4'], ArrayMethods) or
                isinstance(met['p4'], PtEtaPhiMassArrayMethods)):
            raise Exception('Met p4 must be a TLorentzVectorArray!')

    initial_p4 = jet['p4'].copy()  # keep a copy for fixing met
    # initialize the jet momenta to raw values
    _update_jet_ptm(1.0, jet, fromRaw=True)

    # below we work in numpy arrays, JaggedCandidateArray knows how to convert them
    args = {key: getattr(jet, _signature_map[key]).content for key in corrections._jec.signature}
    jec = corrections._jec.getCorrection(**args)

    _update_jet_ptm(jec, jet, fromRaw=True)

    juncs = None
        ## JES
    if corrections._junc is not None:
        args = {key: getattr(jet, _signature_map[key]).content for key in corrections._junc.signature}
        juncs = corrections._junc.getUncertainty(**args)
        juncs = list((name, values) for name, values in juncs)

    # if there's a jer and sf to apply we have to update the momentum too
    # right now only use stochastic smearing
    if corrections._jer is not None and corrections._jersf is not None:
        args = {key: getattr(jet, _signature_map[key]).content for key in corrections._jer.signature}
        jer = corrections._jer.getResolution(**args)

        args = {key: getattr(jet, _signature_map[key]).content for key in corrections._jersf.signature}
        jersf = corrections._jersf.getScaleFactor(**args)

        np.random.seed(10) # sets seed so values from random distributions are reproducible (JER corrections)
        jersmear = jer * np.random.normal(size=jer.size)

        ptGenJet = np.zeros_like(jet.pt.content) if forceStochastic else jet.ptGenJet.content

        doHybrid = ptGenJet > 0

        jsmear = np.zeros(jet._content.size)
        if shift == 'JER_UP':
            jsmear = np.where(doHybrid,
                                 1 + (jersf[:, 1] - 1) * (jet.pt.content - ptGenJet) / jet.pt.content,
                                 1. + np.sqrt(np.maximum(jersf[:, 1]**2 - 1.0, 0)) * jersmear)

        elif shift == 'JER_DW':
            jsmear = np.where(doHybrid,
                                   1 + (jersf[:, -1] - 1) * (jet.pt.content - ptGenJet) / jet.pt.content,
                                   1. + np.sqrt(np.maximum(jersf[:, -1]**2 - 1.0, 0)) * jersmear)
        else:
            jsmear = np.where(doHybrid,
                                  1 + (jersf[:, 0] - 1) * (jet.pt.content - ptGenJet) / jet.pt.content,
                                  1. + np.sqrt(np.maximum(jersf[:, 0]**2 - 1.0, 0)) * jersmear)

        # from PhysicsTools/PatUtils/interface/SmearedJetProducerT.h#L255-L264
        min_jet_pt = corrections.MIN_JET_ENERGY / np.cosh(jet.eta.content)
        min_jet_pt_corr = min_jet_pt / jet.pt.content
        jsmear = np.where(jsmear * jet.pt.content < min_jet_pt,
                              min_jet_pt_corr,
                              jsmear)
        # finally, update the value
        _update_jet_ptm(jsmear, jet)

    # apply JES UP/DW variations
    # have to apply central jersf before calculating junc
    if (corrections._junc is not None) and ((shift == 'JES_UP') or (shift == 'JES_DW')):
        for name, values in juncs:
            jfac = values[:, 0] if shift == 'JES_UP' else values[:, 1]
            _update_jet_ptm(jfac, jet)

    # hack to update the jet p4, we have the fully updated pt and mass here
    jet._content._contents['p4'] = TLorentzVectorArray.from_ptetaphim(jet.pt.content,
                                                                      jet.eta.content,
                                                                      jet.phi.content,
                                                                      jet.mass.content)

    # set MET values
    new_x = met['p4'].x.content - (initial_p4.x - jet['p4'].x).sum()
    new_y = met['p4'].y.content - (initial_p4.y - jet['p4'].y).sum()
    #new_x = new_met_x - (initial_p4.x - jet['p4'].x).sum()
    #new_y = new_met_y - (initial_p4.y - jet['p4'].y).sum()
    _update_met_ptphi(met, pt=np.sqrt(new_x**2 + new_y**2), phi=np.arctan2(new_y, new_x))
    met._content._contents['p4'] = TLorentzVectorArray.from_ptetaphim(
        met.pt.content, met.eta.content, met.phi.content, met.mass.content
    )

    #set_trace()
    if shift == 'MET_UP':
        dx = met._content._contents['MetUnclustEnUpDeltaX']
        dy = met._content._contents['MetUnclustEnUpDeltaY']
        new_met_x = met['p4'].x.content + dx
        new_met_y = met['p4'].y.content + dy
        new_x = new_met_x - (initial_p4.x - jet['p4'].x).sum()
        new_y = new_met_y - (initial_p4.y - jet['p4'].y).sum()
        _update_met_ptphi(met, pt=np.sqrt(new_x**2 + new_y**2), phi=np.arctan2(new_y, new_x))
        met._content._contents['p4'] = TLorentzVectorArray.from_ptetaphim(
            met.pt.content, met.eta.content, met.phi.content, met.mass.content
        )
    elif shift == 'MET_DW':
        dx = met._content._contents['MetUnclustEnUpDeltaX']
        dy = met._content._contents['MetUnclustEnUpDeltaY']
        new_met_x = met['p4'].x.content - dx
        new_met_y = met['p4'].y.content - dy
        new_x = new_met_x - (initial_p4.x - jet['p4'].x).sum()
        new_y = new_met_y - (initial_p4.y - jet['p4'].y).sum()
        _update_met_ptphi(met, pt=np.sqrt(new_x**2 + new_y**2), phi=np.arctan2(new_y, new_x))
        met._content._contents['p4'] = TLorentzVectorArray.from_ptetaphim(
            met.pt.content, met.eta.content, met.phi.content, met.mass.content
        )

    #if corrections._junc is not None:
    #    jets_sin = np.sin(jet['p4'].phi)
    #    jets_cos = np.cos(jet['p4'].phi)
    #    for name, _ in juncs:
    #        for shift in ['up', 'down']:
    #            px = met['p4'].x - (initial_p4.x - jet['pt_{0}_{1}'.format(name, shift)] * jets_cos).sum()
    #            py = met['p4'].y - (initial_p4.y - jet['pt_{0}_{1}'.format(name, shift)] * jets_sin).sum()
    #            met['pt_{0}_{1}'.format(name, shift)] = np.sqrt(px**2 + py**2)
    #            met['phi_{0}_{1}'.format(name, shift)] = np.arctan2(py, px)

