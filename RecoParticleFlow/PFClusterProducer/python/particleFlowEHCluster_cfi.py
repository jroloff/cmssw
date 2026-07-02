import FWCore.ParameterSet.Config as cms

# Default configuration for PFEHClusterProducer.
#
# Input clusters are taken from the standard PF clustering sequence
# (particleFlowClusterECAL / particleFlowClusterHCAL), which means the
# clusters already carry the subdetector-level energy-response calibrations
# applied by PFClusterProducer.  The PF hadronic calibration (combined
# ECAL+HCAL response to pions) is intentionally NOT applied here; it must
# still be applied by PFAlgo when consuming this collection.
#
# Tuning notes
# ------------
# matchingDeltaR = 0.15  is chosen to be slightly larger than a typical
#   HCAL tower size (~0.087) but well below the typical jet radius (0.4),
#   so that we match the HCAL deposit of a single hadron shower without
#   accidentally pulling in deposits from nearby particles.
#
# minEcalEnergy = 0.5 GeV  corresponds to ~12 sigma above the ECAL barrel
#   noise level (~40 MeV/crystal x a few crystals).
#
# minHcalEnergy = 1.0 GeV  corresponds to ~5 sigma above the HCAL noise
#   level (~200 MeV/tower).

pfEHClusterProducer = cms.EDProducer(
    "PFEHClusterProducer",

    # Input PFCluster collections from the standard PF clustering sequence.
    ecalClusters = cms.InputTag("particleFlowClusterECAL"),
    hcalClusters = cms.InputTag("particleFlowClusterHCAL"),

    # Maximum dR in eta-phi between an ECAL and HCAL cluster for them to
    # be merged.  Corresponds roughly to 1.7 HCAL tower widths in the barrel.
    matchingDeltaR = cms.double(0.15),

    # Minimum cluster energies.  Clusters below these thresholds are not used
    # as seeds (ECAL) or candidates for matching (HCAL), though they can
    # still be promoted to unmatched SuperClusters if keepUnmatched=True.
    minEcalEnergy = cms.double(0.5),   # GeV
    minHcalEnergy = cms.double(1.0),   # GeV

    # When True, ECAL or HCAL clusters that find no match in the other
    # subdetector are promoted to single-subdetector SuperClusters so
    # that no energy is dropped from the collection.
    keepUnmatched = cms.bool(True),
)

# -----------------------------------------------------------------------
# Convenience sequence: run after the standard PF cluster producers.
# Insert this before particleFlowBlock if you want to feed EHClusters
# into a custom PFBlockAlgo element type.
# -----------------------------------------------------------------------
from RecoParticleFlow.PFClusterProducer.particleFlowClusterECAL_cff import (
    particleFlowClusterECAL,
)
from RecoParticleFlow.PFClusterProducer.particleFlowClusterHCAL_cff import (
    particleFlowClusterHCAL,
)

pfEHClusterTask = cms.Task(
    particleFlowClusterECAL,
    particleFlowClusterHCAL,
    pfEHClusterProducer,
)

pfEHClusterSequence = cms.Sequence(pfEHClusterTask)

