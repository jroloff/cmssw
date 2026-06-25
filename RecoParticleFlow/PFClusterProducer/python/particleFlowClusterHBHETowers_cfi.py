import FWCore.ParameterSet.Config as cms
from RecoParticleFlow.PFClusterProducer.particleFlowCaloResolution_cfi import _timeResolutionHCALMaxSample
from RecoParticleFlow.PFClusterProducer.particleFlowRecHitHBHEFilters_cfi import *
from RecoParticleFlow.PFClusterProducer.particleFlowClusterHBHE_cfi import *

particleFlowClusterHBHETowers1To26 = particleFlowClusterHBHE.clone(
    recHitsSource = cms.InputTag("particleFlowRecHitHBHEAbsIEta1To26")
)


particleFlowClusterHBHETowers27To29 = particleFlowClusterHBHE.clone(
    recHitsSource = cms.InputTag("particleFlowRecHitHBHEAbsIEta27To29")
)


