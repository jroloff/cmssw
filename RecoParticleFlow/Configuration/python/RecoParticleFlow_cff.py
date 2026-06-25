import FWCore.ParameterSet.Config as cms


from RecoParticleFlow.PFTracking.particleFlowTrack_cff import *
#from RecoParticleFlow.PFTracking.particleFlowTrackWithDisplacedVertex_cff import *

from RecoParticleFlow.PFProducer.particleFlowSimParticle_cff import *
from RecoParticleFlow.PFProducer.particleFlowBlock_cff import *

from RecoParticleFlow.PFProducer.particleFlowEGamma_cff import *
from RecoParticleFlow.PFProducer.particleFlow_cff import *
from RecoParticleFlow.PFProducer.pfElectronTranslator_cff import *
from RecoParticleFlow.PFProducer.pfPhotonTranslator_cff import *
#from RecoParticleFlow.PFProducer.pfGsfElectronCiCSelector_cff import *
from RecoParticleFlow.PFProducer.pfGsfElectronMVASelector_cff import *

from RecoParticleFlow.PFClusterProducer.particleFlowRecHitHBHEFilters_cfi import *
from RecoParticleFlow.PFProducer.pfLinker_cff import *
from RecoParticleFlow.PFClusterProducer.particleFlowClusterHBHETowers_cfi import *
from RecoParticleFlow.PFClusterProducer.particleFlowClusterHCALTowers_cfi import *

from RecoParticleFlow.PFProducer.particleFlowBlockTowers_cfi import *
from RecoParticleFlow.PFProducer.particleFlowTowers_cff import particleFlowTmpTowers

from CommonTools.ParticleFlow.pfParticleSelection_cff import *

from RecoEgamma.EgammaIsolationAlgos.particleBasedIsoProducer_cff import *
from RecoParticleFlow.PFProducer.chargedHadronPFTrackIsolation_cfi import *

from RecoJets.JetProducers.fixedGridRhoProducerFastjet_cfi import *
fixedGridRhoFastjetAllTmp = fixedGridRhoFastjetAll.clone(pfCandidatesTag = "particleFlowTmp")

print("DEBUG: particleFlowTmpTowers =")
particleFlowTmpTask = cms.Task(particleFlowTmp)
particleFlowTmpTowersTask = cms.Task(particleFlowTmpTowers)
particleFlowTmpSeq = cms.Sequence(particleFlowTmpTask)
particleFlowTmpTowersSeq = cms.Sequence(particleFlowTmpTowersTask)

particleFlowRecoTask = cms.Task( particleFlowTrackWithDisplacedVertexTask,
#                                pfGsfElectronCiCSelectionSequence,
                                 pfGsfElectronMVASelectionTask,
                                 particleFlowBlock,
                                 particleFlowEGammaFullTask,
                                 particleFlowTmpTask,
                                 particleFlowTmpTowersTask,
                                 fixedGridRhoFastjetAllTmp,
                                 particleFlowTmpPtrs,
                                 particleFlowEGammaFinalTask,
                                 pfParticleSelectionTask )
particleFlowReco = cms.Sequence(particleFlowRecoTask)

particleFlowLinksTask = cms.Task( particleFlow, particleFlowPtrs, chargedHadronPFTrackIsolation, particleBasedIsolationTask)
particleFlowLinks = cms.Sequence(particleFlowLinksTask)

# In your custom cff or via a modifier:
particleFlowTowers_Task = cms.Task(
    particleFlowRecHitHBHEAbsIEta1To26,
    particleFlowRecHitHBHEAbsIEta27To29,
    particleFlowClusterHBHETowers1To26,
    particleFlowClusterHBHETowers27To29,
    particleFlowClusterHCALTowers,
    particleFlowClusterHCALTowers1To26,
    particleFlowClusterHCALTowers27To29,
    particleFlowBlockTowers,
    particleFlowTmpTowers,
    #particleFlowTowers      # your final PFCandidate collection
)
particleFlowRecoTowers = cms.Sequence(particleFlowTowers_Task)

#
# for phase 2
particleFlowTmpBarrel = particleFlowTmp.clone()
_phase2_hgcal_particleFlowTmp = cms.EDProducer(
    "PFCandidateListMerger",
    src = cms.VInputTag("particleFlowTmpBarrel",
                        "pfTICL")

)

from Configuration.Eras.Modifier_phase2_hgcal_cff import phase2_hgcal
phase2_hgcal.toReplaceWith( particleFlowTmp, _phase2_hgcal_particleFlowTmp )
phase2_hgcal.toModify(
    particleFlowTmpBarrel,
    vetoEndcap = True
    # If true, PF(Muon)Algo will ignore muon candidates incorporated via pfTICL
    # in addMissingMuons. This will prevent potential double-counting.
)

#
# for simPF
from RecoParticleFlow.PFClusterProducer.particleFlowClusterHGC_cfi import *
from RecoParticleFlow.PFTracking.hgcalTrackCollection_cfi import *
from RecoParticleFlow.PFProducer.simPFProducer_cff import *
from SimTracker.TrackerHitAssociation.tpClusterProducer_cfi import *
from SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi import *

_phase2_hgcal_simPFTask = cms.Task( pfTrack ,
                                    hgcalTrackCollection ,
                                    tpClusterProducer ,
                                    quickTrackAssociatorByHits ,
                                    particleFlowClusterHGCalFromSimCl ,
                                    simPFProducer )
_phase2_hgcal_simPFSequence = cms.Sequence(_phase2_hgcal_simPFTask)

_phase2_hgcal_particleFlowRecoTask = cms.Task( _phase2_hgcal_simPFTask , particleFlowRecoTask.copy() )
_phase2_hgcal_particleFlowRecoTask.replace( particleFlowTmpTask, cms.Task( particleFlowTmpBarrel, particleFlowTmp ) )

phase2_hgcal.toReplaceWith( particleFlowRecoTask, _phase2_hgcal_particleFlowRecoTask )

#
# for heavy ion
from Configuration.Eras.Modifier_pp_on_XeXe_2017_cff import pp_on_XeXe_2017
from Configuration.ProcessModifiers.pp_on_AA_cff import pp_on_AA

for e in [pp_on_XeXe_2017, pp_on_AA]:
    e.toModify(particleFlowDisplacedVertexCandidate,
               tracksSelectorParameters = dict(pt_min = 999999.0,
                                               nChi2_max = 0.0,
                                               pt_min_prim = 999999.0,
                                               dxy = 999999.0)
               )
    e.toModify(pfNoPileUpIso, enable = False)
    e.toModify(pfPileUpIso, enable = False)
    e.toModify(pfNoPileUp, enable = False)
    e.toModify(pfPileUp, enable = False)


# for MLPF, replace standard PFAlgo with the ONNX-based MLPF producer 
from Configuration.ProcessModifiers.mlpf_cff import mlpf
from RecoParticleFlow.PFProducer.mlpfProducer_cfi import mlpfProducer
mlpf.toReplaceWith(particleFlowTmp, mlpfProducer)

#
# switch from pfTICL to simPF
def _findIndicesByModule(process,name):
    ret = []
    if hasattr(process,'particleFlowBlock'):
        for i, pset in enumerate(process.particleFlowBlock.elementImporters):
            if pset.importerName.value() == name:
                ret.append(i)
    return ret

from Configuration.Eras.Modifier_phase2_hgcal_cff import phase2_hgcal
def replaceTICLwithSimPF(process):
    if hasattr(process,'particleFlowTmp'):
      process.particleFlowTmp.src = ['particleFlowTmpBarrel', 'simPFProducer']

    if hasattr(process,'particleFlowTmpBarrel'):
      process.particleFlowTmpBarrel.vetoEndcap = False

    _insertTrackImportersWithVeto = {}
    _trackImporters = ['GeneralTracksImporter','ConvBremTrackImporter',
                   'ConversionTrackImporter','NuclearInteractionTrackImporter']
    for importer in _trackImporters:
        for idx in _findIndicesByModule(process,importer):
            _insertTrackImportersWithVeto[idx] = dict(
                vetoMode = cms.uint32(0), # HGCal-region PFTrack list for simPF
                vetoSrc = cms.InputTag('hgcalTrackCollection:TracksInHGCal')
            )
    phase2_hgcal.toModify(
      process.particleFlowBlock,
      elementImporters = _insertTrackImportersWithVeto
    )

    return process
