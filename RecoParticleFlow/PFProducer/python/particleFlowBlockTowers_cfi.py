from RecoParticleFlow.PFProducer.particleFlowBlock_cfi import particleFlowBlock
import FWCore.ParameterSet.Config as cms

particleFlowBlockTowers = particleFlowBlock.clone()

for i, pset in enumerate(particleFlowBlockTowers.elementImporters):
    if pset.importerName.value() == 'GenericClusterImporter':
        if pset.source.value() == 'particleFlowClusterHCAL':
            particleFlowBlockTowers.elementImporters[i].source = cms.InputTag(
                "particleFlowClusterHCALTowers"
            )

