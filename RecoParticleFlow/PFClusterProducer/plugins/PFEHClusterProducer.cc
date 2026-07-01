// PFEHClusterProducer
// --------------------------------
// EDProducer that merges reco::PFCluster objects from the ECAL and HCAL into
// combined PFEHCluster objects by associating clusters that are
// within a configurable cone (deltaR) in eta-phi space.
//
// Algorithm
// ---------
// 1. For every ECAL cluster, find all HCAL clusters within dR < matchingDeltaR_.
//    Multiple HCAL clusters can be associated to one ECAL seed (sub-clusters
//    within the cone are all merged), and one HCAL cluster can be shared by
//    at most one ECAL seed (it is assigned to the nearest ECAL seed).
// 2. ECAL clusters with no HCAL match and HCAL clusters with no ECAL match
//    are promoted to an EHCluster 
// 3. Position is recomputed as the energy-weighted centroid of all
//    constituent cluster positions.
//
// Sharing / ambiguity resolution
// --------------------------------
// When two ECAL clusters both fall within dR of the same HCAL cluster the
// HCAL cluster is assigned exclusively to the nearest ECAL seed.  This is
// consistent with the philosophy used in PFBlockAlgo for track-cluster links.

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFEHCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <vector>

class PFEHClusterProducer : public edm::stream::EDProducer<> {
public:
  explicit PFEHClusterProducer(const edm::ParameterSet& iConfig);
  ~PFEHClusterProducer() override = default;

  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // ------------------------------------------------------------------ helpers

  // Recompute the energy-weighted centroid from all constituent clusters.
  static math::XYZPoint energyWeightedPosition(
      const reco::PFClusterRefVector& ecalClusters,
      const reco::PFClusterRefVector& hcalClusters);

  const edm::EDPutTokenT<reco::PFEHClusterCollection> pfehClusterToken_;
  // ------------------------------------------------------------------ tokens
  const edm::EDGetTokenT<reco::PFClusterCollection> ecalToken_;
  const edm::EDGetTokenT<reco::PFClusterCollection> hcalToken_;


  // ------------------------------------------------------------------ config
  const double matchingDeltaR_;   // cone half-opening angle for ECAL-HCAL matching
};

// ============================================================================
// Constructor
// ============================================================================
PFEHClusterProducer::PFEHClusterProducer(
    const edm::ParameterSet& iConfig)
    : pfehClusterToken_{produces<reco::PFEHClusterCollection>()}, 
      ecalToken_(consumes<reco::PFClusterCollection>(
          iConfig.getParameter<edm::InputTag>("ecalClusters"))),
      hcalToken_(consumes<reco::PFClusterCollection>(
          iConfig.getParameter<edm::InputTag>("hcalClusters"))),
      matchingDeltaR_(iConfig.getParameter<double>("matchingDeltaR")){
}

// ============================================================================
// produce
// ============================================================================
void PFEHClusterProducer::produce(edm::Event& iEvent,
                                              const edm::EventSetup& /*iSetup*/) {
  // --- retrieve input collections -------------------------------------------
  edm::Handle<reco::PFClusterCollection> ecalHandle;
  iEvent.getByToken(ecalToken_, ecalHandle);
  // TODO need some sort of error handling

  edm::Handle<reco::PFClusterCollection> hcalHandle;
  iEvent.getByToken(hcalToken_, hcalHandle);

  const auto& ecalClusters = *ecalHandle;
  const auto& hcalClusters = *hcalHandle;

  const int nEcal = ecalClusters.size();
  const int nHcal = hcalClusters.size();

  // --- output collection ----------------------------------------------------
  auto output = std::make_unique<reco::PFEHClusterCollection>();
  output->reserve(nEcal);  // at most one SC per ECAL cluster

  // --- bookkeeping: track which HCAL clusters have been assigned ------------
  // hcalAssignedTo[j] = index of the ECAL seed cluster that owns HCAL cluster j
  // -1 means unassigned.
  std::vector<int> hcalAssignedTo(nHcal, -1);

  // For each HCAL cluster, find the nearest ECAL seed within matchingDeltaR_.
  // This nearest-neighbour pass resolves sharing unambiguously.
  std::vector<double> hcalMinDR2(nHcal, std::numeric_limits<double>::max());

  const double dR2cut = matchingDeltaR_ * matchingDeltaR_;

  // First pass: determine ownership of each HCAL cluster.
  for (int iEcal = 0; iEcal < nEcal; ++iEcal) {
    const auto& ecalClus = ecalClusters[iEcal];

    const double ecalEta = ecalClus.eta();
    const double ecalPhi = ecalClus.phi();

    for (int iHcal = 0; iHcal < nHcal; ++iHcal) {
      const auto& hcalClus = hcalClusters[iHcal];

      const double dR2 = reco::deltaR2(ecalEta, ecalPhi,
                                        hcalClus.eta(), hcalClus.phi());
      if (dR2 < dR2cut && dR2 < hcalMinDR2[iHcal]) {
        hcalMinDR2[iHcal]    = dR2;
        hcalAssignedTo[iHcal] = iEcal;
      }
    }
  }

  // Second pass: build one EHCluster per ECAL seed.
  // Track which ECAL clusters were actually used (for unmatched handling).
  std::vector<bool> ecalUsed(nEcal, false);

  for (int iEcal = 0; iEcal < nEcal; ++iEcal) {
    const auto& ecalClus = ecalClusters[iEcal];

    // Collect all HCAL clusters assigned to this ECAL seed.
    reco::PFClusterRefVector hcalRefs;
    for (int iHcal = 0; iHcal < nHcal; ++iHcal) {
      if (hcalAssignedTo[iHcal] == iEcal) {
        hcalRefs.push_back(reco::PFClusterRef(hcalHandle, iHcal));
      }
    }

    // Collect the ECAL ref for this seed.
    reco::PFClusterRefVector ecalRefs;
    reco::PFClusterRef seedRef(ecalHandle, iEcal);
    ecalRefs.push_back(seedRef);

    // Compute energies.
    double eEcal = ecalClus.energy();
    double eHcal = 0.;
    for (const auto& href : hcalRefs) eHcal += href->energy();

    // Build the EHCluster.
    reco::PFEHCluster ehClust(seedRef, ecalRefs, hcalRefs);
    ehClust.setRawEcalEnergy(eEcal);
    ehClust.setRawHcalEnergy(eHcal);
    ehClust.setPosition(energyWeightedPosition(ecalRefs, hcalRefs));

    output->push_back(std::move(ehClust));
    ecalUsed[iEcal] = true;

    LogDebug("PFEHClusterProducer")
        << "SC seed ecalE=" << eEcal << " GeV, nHcal=" << hcalRefs.size()
        << " hcalE=" << eHcal << " GeV"
        << " eta=" << output->back().eta() << " phi=" << output->back().phi();
  }

  // Promote unmatched HCAL clusters to HCAL-only SCs.
  for (int iHcal = 0; iHcal < nHcal; ++iHcal) {
    if (hcalAssignedTo[iHcal] >= 0) continue;  // already assigned
    const auto& hcalClus = hcalClusters[iHcal];

    reco::PFClusterRef hRef(hcalHandle, iHcal);
    reco::PFClusterRefVector emptyEcal;
    reco::PFClusterRefVector hcalRefs;
    hcalRefs.push_back(hRef);

    reco::PFEHCluster ehClust(hRef, emptyEcal, hcalRefs);
    ehClust.setRawEcalEnergy(0.);
    ehClust.setRawHcalEnergy(hcalClus.energy());
    ehClust.setPosition(math::XYZPoint(hcalClus.position().x(),
                                  hcalClus.position().y(),
                                  hcalClus.position().z()));
    output->push_back(std::move(ehClust));
  }

  // Keep unmatched ECAL-only clusters.
  for (int iEcal = 0; iEcal < nEcal; ++iEcal) {
    if (ecalUsed[iEcal]) continue;
    const auto& ecalClus = ecalClusters[iEcal];

    reco::PFClusterRef eRef(ecalHandle, iEcal);
    reco::PFClusterRefVector ecalRefs;
    ecalRefs.push_back(eRef);
    reco::PFClusterRefVector emptyHcal;

    reco::PFEHCluster ehClust(eRef, ecalRefs, emptyHcal);
    ehClust.setRawEcalEnergy(ecalClus.energy());
    ehClust.setRawHcalEnergy(0.);
    ehClust.setPosition(math::XYZPoint(ecalClus.position().x(),
                                  ecalClus.position().y(),
                                  ecalClus.position().z()));
    output->push_back(std::move(ehClust));
  }
  std::cout << output->size() << "\t" << nEcal << "\t" << nHcal << std::endl;

  LogDebug("PFEHClusterProducer")
      << "Produced " << output->size() << " PFEHClusters from "
      << nEcal << " ECAL and " << nHcal << " HCAL PFClusters.";

  iEvent.put(pfehClusterToken_, std::move(output));
}

// ============================================================================
// energyWeightedPosition
// ============================================================================
math::XYZPoint PFEHClusterProducer::energyWeightedPosition(
    const reco::PFClusterRefVector& ecalClusters,
    const reco::PFClusterRefVector& hcalClusters) {
  double sumW = 0., sumX = 0., sumY = 0., sumZ = 0.;

  auto accumulate = [&](const reco::PFClusterRefVector& clusters) {
    for (const auto& ref : clusters) {
      const double w = ref->energy();
      const auto& pos = ref->position();
      sumW += w;
      sumX += w * pos.x();
      sumY += w * pos.y();
      sumZ += w * pos.z();
    }
  };

  accumulate(ecalClusters);
  accumulate(hcalClusters);

  if (sumW <= 0.) return math::XYZPoint(0., 0., 0.);
  return math::XYZPoint(sumX / sumW, sumY / sumW, sumZ / sumW);
}

// ============================================================================
// fillDescriptions
// ============================================================================
void PFEHClusterProducer::fillDescriptions(
    edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("ecalClusters", edm::InputTag("particleFlowClusterECAL"))
      ->setComment("ECAL PFCluster collection");
  desc.add<edm::InputTag>("hcalClusters", edm::InputTag("particleFlowClusterHCAL"))
      ->setComment("HCAL PFCluster collection");
  desc.add<double>("matchingDeltaR", 0.15)
      ->setComment("Maximum dR between ECAL and HCAL cluster positions for matching");
  descriptions.add("pfEHClusterProducer", desc);
}

// ============================================================================
// Plugin registration
// ============================================================================
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFEHClusterProducer);
