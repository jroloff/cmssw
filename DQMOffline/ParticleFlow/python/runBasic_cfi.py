import FWCore.ParameterSet.Config as cms
from DQMOffline.ParticleFlow.basicConfig_cff import *
from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer

PFAnalyzer = DQMEDAnalyzer("PFAnalyzer",
    pfJetCollection        = cms.InputTag("ak4PFJetsCHS"),
    pfCandidates             = cms.InputTag("particleFlow"),
    PVCollection             = cms.InputTag("offlinePrimaryVerticesWithBS"),


    jetAnalysis2 = jetAnalysis.clone(),


)
