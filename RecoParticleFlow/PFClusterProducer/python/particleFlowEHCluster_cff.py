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
#from RecoParticleFlow.PFClusterProducer.pfEHClusterProducer_cfi import (
#    pfEHClusterProducer,
#)

pfEHClusterProducer = cms.EDProducer(
    "PFEHClusterProducer",
    ## Input PFCluster collections from the standard PF clustering sequence.
    #ecalClusters = cms.InputTag("particleFlowClusterECAL"),
    #hcalClusters = cms.InputTag("particleFlowClusterHCAL"),
    ## Maximum dR in eta-phi between an ECAL and HCAL cluster for them to
    ## be merged.  Corresponds roughly to 1.7 HCAL tower widths in the barrel.
    #matchingDeltaR = cms.double(0.2),
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
#from RecoParticleFlow.PFProducer.particleFlowBlock_cfi import particleFlowBlock
#
#pfEHBlock = particleFlowBlock.clone(
#    # Replace the bare ECAL and HCAL cluster sources with the merged EH
#    # collection.  The elementImporters list controls which collections
#    # PFBlockAlgo reads; we keep everything except the bare ECAL/HCAL
#    # cluster importers and substitute a new EHCluster importer.
#    #
#    # In the standard particleFlowBlock the importers are:
#    #   [GSFTrackImporter, ConvBremTrackImporter, GeneralTracksImporter,
#    #    SuperClusterImporter, ECALClusterImporter, HCALClusterImporter, ...]
#    #
#    # We add an EHClusterImporter entry.  The existing ECAL and HCAL
#    # importers are kept so that the standard PF algorithm also runs
#    # unmodified in the same block job; downstream, PFEHCandidateProducer
#    # selects only the EH-type elements.
#    elementImporters = particleFlowBlock.elementImporters + [
#        cms.PSet(
#            importerName = cms.string("EHClusterImporter"),
#            source       = cms.InputTag("pfEHClusterProducer"),
#            # BES (block element source) key must match the string
#            # registered in the importer plugin.
#            importerDef  = cms.PSet(),
#        )
#    ],
#
#    # Append the new Track-EHCluster linker to the existing linker list.
#    # PFBlockAlgo iterates over this list and calls each linker for every
#    # pair of element types it handles.
#    # Append the KDTree-based Track-EHCluster linker.
#    # The linkerName must match DEFINE_EDM_PLUGIN in KDTreeLinkerTrackEH.cc.
#    # Parameters mirror the defaults of KDTreeLinkerTrackEcal (nSigmaECAL)
#    # and KDTreeLinkerTrackHcal (nSigmaHCAL, crystalSize*).
#    linkDefinitions = particleFlowBlock.linkDefinitions + [
#        cms.PSet(
#            linkerName     = cms.string("KDTreeLinkerTrackAndEH"),
#            linkType       = cms.string("TRACK:EH"),
#            useKDTree      = cms.bool(True),
#            nSigmaECAL     = cms.double(2.0),
#            nSigmaHCAL     = cms.double(1.0),
#            crystalSizeEta = cms.double(0.04),
#            crystalSizePhi = cms.double(0.04),
#        )
#    ],
#)

# -----------------------------------------------------------------------
# Step 3: create PFCandidates from the blocks.
# Point at the cloned block producer, not the standard one.
# -----------------------------------------------------------------------
#particleFlowEH = cms.EDProducer(
#    "PFEHCandidateProducer",
#    pfBlocks               = cms.InputTag("pfEHBlock"),
#    calibAlpha             = cms.double(1.00),
#    calibBeta              = cms.double(1.30),
#    nSigmaExcess           = cms.double(1.0),
#    minPhotonEnergy        = cms.double(0.5),
#    minNeutralHadronEnergy = cms.double(1.0),
#)

# -----------------------------------------------------------------------
# Task and Sequence
# -----------------------------------------------------------------------

particleFlowEHTask = cms.Task(
    pfEHClusterProducer,
    #pfEHBlock,
    #particleFlowEH,
)

particleFlowEHCluster = cms.Sequence(particleFlowEHTask)



