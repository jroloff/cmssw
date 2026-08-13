import FWCore.ParameterSet.Config as cms
import copy

# particleFlowEHCluster_cff
# -------------------------
# Defines the full task chain that produces a reco::PFCandidateCollection
# built from PFEHClusters, using the standard PFBlockAlgo machinery with
# a new KDTreeLinkerTrackAndEH plugged in.
#
# Sequence:
#
#   particleFlowClusterECAL  \
#                             +--> pfEHClusterProducer
#   particleFlowClusterHCAL  /         |
#                                       v
#                              pfEHBlock   (standard PFBlockAlgo, with
#                                           KDTreeLinkerTrackAndEH added)
#                                       |
#                                       v
#                              particleFlowEH   (PFEHCandidateProducer)
#
# The output label "particleFlowEH" is distinct from "particleFlow" so both
# collections can coexist in the same event for validation.

# -----------------------------------------------------------------------
# Step 1: merge ECAL + HCAL PFClusters into PFEHClusters
# -----------------------------------------------------------------------
pfEHClusterProducer = cms.EDProducer(
    "PFEHClusterProducer",
    ## Input PFCluster collections from the standard PF clustering sequence.
    #ecalClusters = cms.InputTag("particleFlowClusterECAL"),
    #hcalClusters = cms.InputTag("particleFlowClusterHCAL"),
    ## Maximum dR in eta-phi between an ECAL and HCAL cluster for them to
    ## be merged.  Corresponds roughly to 1.7 HCAL tower widths in the barrel.
    #matchingDeltaR = cms.double(0.0),
    # Input PFCluster collections from the standard PF clustering sequence.
    ecalClusters = cms.InputTag("particleFlowClusterECAL"),
    hcalClusters = cms.InputTag("particleFlowClusterHCAL"),
    hfClusters   = cms.InputTag("particleFlowClusterHF"),
    hoClusters   = cms.InputTag("particleFlowClusterHO"),
    # Maximum dR in eta-phi between an ECAL cluster and an HCAL/HF/HO
    # cluster for them to be merged into a single PFEHCluster.
    matchingDeltaR   = cms.double(0.0),
    matchingDeltaRHF = cms.double(0.0),
    matchingDeltaRHO = cms.double(0.0),
    #matchingDeltaR   = cms.double(0.15),
    #matchingDeltaRHF = cms.double(0.3),
    #matchingDeltaRHO = cms.double(0.2),
)

# -----------------------------------------------------------------------
# Step 2: run PFBlockAlgo with the EH linker added.
#
# We clone the standard particleFlowBlock producer and:
#   a) replace its element sources so it sees PFEHClusters instead of (or
#      in addition to) bare ECAL/HCAL PFClusters, and
#   b) append KDTreeLinkerTrackAndEH to its list of linkers.
#
# The clone means we inherit all the existing TRACK-ECAL, TRACK-HCAL,
# ECAL-HCAL, muon, GSF, etc. linkers for free — we only add one new one.
# -----------------------------------------------------------------------
from RecoParticleFlow.PFProducer.particleFlowBlock_cfi import particleFlowBlock

_excludedImporters = {"ECALClusterImporter", "HCALClusterImporter", "HOClusterImporter", "HFClusterImporter"}

_replacedClusterImporters = {
    ("ECALClusterImporter", "particleFlowClusterECAL"),
    ("GenericClusterImporter", "particleFlowClusterHCAL"),
    ("GenericClusterImporter", "particleFlowClusterHF"),
    ("GenericClusterImporter", "particleFlowClusterHO"),
}
 
 
def _isReplacedImporter(imp):
    if not hasattr(imp, "source"):
        return False
    return (imp.importerName.value(), imp.source.getModuleLabel()) in _replacedClusterImporters
 
 
pfEHBlock = particleFlowBlock.clone(
    elementImporters = [imp for imp in particleFlowBlock.elementImporters if not _isReplacedImporter(imp)] + [
        # Re-add the ECAL importer, now sourced from the unmerged-ECAL
        # filtered collection. NOTE: ECALClusterImporter also needs a
        # BCtoPFCMap (SuperCluster <-> PFCluster association ValueMap) that
        # is keyed off the *original* particleFlowClusterECAL collection.
        # If SuperCluster matching against unmerged-only ECAL clusters turns
        # out to behave differently than expected, that ValueMap is the
        # first place to look.
        cms.PSet(
            importerName = cms.string("ECALClusterImporter"),
            source       = cms.InputTag("pfEHClusterProducer", "unmergedECAL"),
            BCtoPFCMap   = cms.InputTag("particleFlowSuperClusterECAL", "PFClusterAssociationEBEE"),
        ),
        cms.PSet(
            importerName = cms.string("GenericClusterImporter"),
            source       = cms.InputTag("pfEHClusterProducer", "unmergedHCAL"),
        ),
        cms.PSet(
            importerName = cms.string("GenericClusterImporter"),
            source       = cms.InputTag("pfEHClusterProducer", "unmergedHF"),
        ),
        cms.PSet(
            importerName = cms.string("GenericClusterImporter"),
            source       = cms.InputTag("pfEHClusterProducer", "unmergedHO"),
        ),
        # Merged EH clusters.
        cms.PSet(
            importerName = cms.string("EHClusterImporter"),
            source       = cms.InputTag("pfEHClusterProducer"),
        ),
    ],
 
    # Append the new Track-EHCluster linker to the existing linker list.
    # PFBlockAlgo iterates over this list and calls each linker for every
    # pair of element types it handles.
    #
    # linkerName must match the DEFINE_EDM_PLUGIN registration in
    # TrackAndEHClusterLinker.cc, i.e. "TrackAndEHClusterLinker" (NOT
    # "KDTreeLinkerTrackAndEH"). When useKDTree=True, PFBlockAlgo looks up
    # the companion KDTree linker under the name "KDTree" + linkerName, i.e.
    # "KDTreeTrackAndEHClusterLinker" -- which matches the registration in
    # KDTreeLinkerTrackEHCluster.cc. TrackAndEHClusterLinker itself only
    # reads "useKDTree" (and an untracked "debug"); it hardcodes the HCAL
    # entrance/exit trajectory points internally, so no
    # trajectoryLayerEntrance/trajectoryLayerExit parameters are needed here.
    linkDefinitions = particleFlowBlock.linkDefinitions + [
        cms.PSet(
            linkerName = cms.string("TrackAndEHClusterLinker"),
            linkType   = cms.string("TRACK:EH"),
            useKDTree  = cms.bool(True),
            trajectoryLayerEntrance = cms.string("HCALEntrance"),
            trajectoryLayerExit     = cms.string("HCALExit"),
        ),
    ],
)



# -----------------------------------------------------------------------
# Step 3: create PFCandidates from the blocks.
# Point at the cloned block producer, not the standard one.
# -----------------------------------------------------------------------
#from RecoParticleFlow.PFProducer.particleFlow_cff import *
from RecoParticleFlow.PFProducer.particleFlow_cfi import particleFlow as _particleFlow

particleFlowEH = _particleFlow.clone(blocks = cms.InputTag("pfEHBlock"))

# -----------------------------------------------------------------------
# Task and Sequence
# -----------------------------------------------------------------------

particleFlowEHTask = cms.Task(
    pfEHClusterProducer,
    pfEHBlock,
    particleFlowEH,
)

particleFlowEHCluster = cms.Sequence(particleFlowEHTask)

