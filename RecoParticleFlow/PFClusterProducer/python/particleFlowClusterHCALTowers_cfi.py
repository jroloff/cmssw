import FWCore.ParameterSet.Config as cms
from RecoParticleFlow.PFClusterProducer.particleFlowClusterHCAL_cfi import particleFlowClusterHCAL

particleFlowClusterHCALTowers1To26 = particleFlowClusterHCAL.clone(
    clustersSource = cms.InputTag("particleFlowClusterHBHETowers1To26")
)

particleFlowClusterHCALTowers27To29 = particleFlowClusterHCAL.clone(
    clustersSource = cms.InputTag("particleFlowClusterHBHETowers27To29")
)

particleFlowClusterHCALTowers = cms.EDProducer(
    "PFClusterCollectionMerger",
    inputs = cms.VInputTag(
        cms.InputTag("particleFlowClusterHCALTowers1To26"),
        cms.InputTag("particleFlowClusterHCALTowers27To29"),
    )
)



