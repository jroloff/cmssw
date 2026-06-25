from RecoParticleFlow.PFProducer.particleFlow_cff import particleFlowTmp
import FWCore.ParameterSet.Config as cms


particleFlowTmpTowers = particleFlowTmp.clone(
    blocks = cms.InputTag("particleFlowBlockTowers")
)

