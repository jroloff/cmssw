// PFEHClusterProducer
// --------------------------------
// EDProducer that merges reco::PFCluster objects from the ECAL, HCAL, HF, and
// HO into combined PFEHCluster objects by associating clusters that are
// within a configurable cone (deltaR) in eta-phi space.
//
// Algorithm
// ---------
// 1. For every ECAL cluster, find all HCAL/HF/HO clusters within dR < the
//    matching cone for that subdetector. Multiple clusters from a given
//    subdetector can be associated to one ECAL seed (sub-clusters within the
//    cone are all merged), and one HCAL/HF/HO cluster can be shared by at
//    most one ECAL seed (it is assigned to the nearest ECAL seed).
// 2. ECAL clusters with no match in any subdetector, and HCAL/HF/HO clusters
//    with no ECAL match, are promoted to their own EHCluster.
// 3. Position is recomputed as the energy-weighted centroid of all
//    constituent cluster positions.
//
// Sharing / ambiguity resolution
// --------------------------------
// When two ECAL clusters both fall within dR of the same HCAL/HF/HO cluster,
// that cluster is assigned exclusively to the nearest ECAL seed. This is
// consistent with the philosophy used in PFBlockAlgo for track-cluster links,
// and is applied independently per subdetector (an HCAL cluster's assignment
// does not affect HF or HO assignment).

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

#include <limits>
#include <vector>

class PFEHClusterProducer : public edm::stream::EDProducer<> {
public:
  explicit PFEHClusterProducer(const edm::ParameterSet& iConfig);
  ~PFEHClusterProducer() override = default;

  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // ------------------------------------------------------------------ helpers

  // For every cluster in `candidates`, find the nearest ECAL seed within
  // `dR2cut` and record ownership in `assignedTo` (index into ecalClusters,
  // -1 if no match). `minDR2` is working storage reused across calls.
  static void assignNearestEcalSeed(const reco::PFClusterCollection& ecalClusters,
                                     const reco::PFClusterCollection& candidates,
                                     double dR2cut,
                                     std::vector<int>& assignedTo,
                                     std::vector<double>& minDR2);

  // Recompute the energy-weighted centroid from all constituent clusters.
  static math::XYZPoint energyWeightedPosition(
      const reco::PFClusterRefVector& ecalClusters,
      const reco::PFClusterRefVector& hcalClusters,
      const reco::PFClusterRefVector& hfClusters,
      const reco::PFClusterRefVector& hoClusters);

  const edm::EDPutTokenT<reco::PFEHClusterCollection> pfehClusterToken_;
  const edm::EDPutTokenT<reco::PFClusterCollection> unmergedEcalToken_;
  const edm::EDPutTokenT<reco::PFClusterCollection> unmergedHcalToken_;
  const edm::EDPutTokenT<reco::PFClusterCollection> unmergedHfToken_;
  const edm::EDPutTokenT<reco::PFClusterCollection> unmergedHoToken_;

  // ------------------------------------------------------------------ tokens
  const edm::EDGetTokenT<reco::PFClusterCollection> ecalToken_;
  const edm::EDGetTokenT<reco::PFClusterCollection> hcalToken_;
  const edm::EDGetTokenT<reco::PFClusterCollection> hfToken_;
  const edm::EDGetTokenT<reco::PFClusterCollection> hoToken_;

  // ------------------------------------------------------------------ config
  const double matchingDeltaR_;    // cone half-opening angle for ECAL-HCAL matching
  const double matchingDeltaRHF_;  // cone half-opening angle for ECAL-HF matching
  const double matchingDeltaRHO_;  // cone half-opening angle for ECAL-HO matching
};

// ============================================================================
// Constructor
// ============================================================================
PFEHClusterProducer::PFEHClusterProducer(
    const edm::ParameterSet& iConfig)
    : pfehClusterToken_{produces<reco::PFEHClusterCollection>()},
      unmergedEcalToken_{produces<reco::PFClusterCollection>("unmergedECAL")},
      unmergedHcalToken_{produces<reco::PFClusterCollection>("unmergedHCAL")},
      unmergedHfToken_{produces<reco::PFClusterCollection>("unmergedHF")},
      unmergedHoToken_{produces<reco::PFClusterCollection>("unmergedHO")},
      ecalToken_(consumes<reco::PFClusterCollection>(
          iConfig.getParameter<edm::InputTag>("ecalClusters"))),
      hcalToken_(consumes<reco::PFClusterCollection>(
          iConfig.getParameter<edm::InputTag>("hcalClusters"))),
      hfToken_(consumes<reco::PFClusterCollection>(
          iConfig.getParameter<edm::InputTag>("hfClusters"))),
      hoToken_(consumes<reco::PFClusterCollection>(
          iConfig.getParameter<edm::InputTag>("hoClusters"))),
      matchingDeltaR_(iConfig.getParameter<double>("matchingDeltaR")),
      matchingDeltaRHF_(iConfig.getParameter<double>("matchingDeltaRHF")),
      matchingDeltaRHO_(iConfig.getParameter<double>("matchingDeltaRHO")) {
}

// ============================================================================
// assignNearestEcalSeed
// ============================================================================
void PFEHClusterProducer::assignNearestEcalSeed(
    const reco::PFClusterCollection& ecalClusters,
    const reco::PFClusterCollection& candidates,
    double dR2cut,
    std::vector<int>& assignedTo,
    std::vector<double>& minDR2) {
  const int nEcal = ecalClusters.size();
  const int nCand = candidates.size();

  assignedTo.assign(nCand, -1);
  minDR2.assign(nCand, std::numeric_limits<double>::max());

  for (int iEcal = 0; iEcal < nEcal; ++iEcal) {
    const auto& ecalClus = ecalClusters[iEcal];
    const double ecalEta = ecalClus.eta();
    const double ecalPhi = ecalClus.phi();

    for (int iCand = 0; iCand < nCand; ++iCand) {
      const auto& candClus = candidates[iCand];

      const double dR2 = reco::deltaR2(ecalEta, ecalPhi,
                                        candClus.eta(), candClus.phi());
      if (dR2 < dR2cut && dR2 < minDR2[iCand]) {
        minDR2[iCand] = dR2;
        assignedTo[iCand] = iEcal;
      }
    }
  }
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

  edm::Handle<reco::PFClusterCollection> hfHandle;
  iEvent.getByToken(hfToken_, hfHandle);

  edm::Handle<reco::PFClusterCollection> hoHandle;
  iEvent.getByToken(hoToken_, hoHandle);

  const auto& ecalClusters = *ecalHandle;
  const auto& hcalClusters = *hcalHandle;
  const auto& hfClusters = *hfHandle;
  const auto& hoClusters = *hoHandle;

  const int nEcal = ecalClusters.size();
  const int nHcal = hcalClusters.size();
  const int nHf = hfClusters.size();
  const int nHo = hoClusters.size();

  // --- output collection ----------------------------------------------------
  auto output = std::make_unique<reco::PFEHClusterCollection>();
  output->reserve(nEcal+nHcal + nHf + nHo); // At most, one cluster per cluster


  auto unmergedEcalOut = std::make_unique<reco::PFClusterCollection>();
  auto unmergedHcalOut = std::make_unique<reco::PFClusterCollection>();
  auto unmergedHfOut = std::make_unique<reco::PFClusterCollection>();
  auto unmergedHoOut = std::make_unique<reco::PFClusterCollection>();
 

  // --- ownership passes: each candidate collection assigned independently --
  std::vector<int> hcalAssignedTo, hfAssignedTo, hoAssignedTo;
  std::vector<double> scratchDR2;  // reused working buffer

  assignNearestEcalSeed(ecalClusters, hcalClusters,
                         matchingDeltaR_ * matchingDeltaR_, hcalAssignedTo, scratchDR2);
  assignNearestEcalSeed(ecalClusters, hfClusters,
                         matchingDeltaRHF_ * matchingDeltaRHF_, hfAssignedTo, scratchDR2);
  assignNearestEcalSeed(ecalClusters, hoClusters,
                         matchingDeltaRHO_ * matchingDeltaRHO_, hoAssignedTo, scratchDR2);

  // Second pass: build one EHCluster per ECAL seed.
  std::vector<bool> ecalUsed(nEcal, false);

  for (int iEcal = 0; iEcal < nEcal; ++iEcal) {
    const auto& ecalClus = ecalClusters[iEcal];

    reco::PFClusterRefVector hcalRefs, hfRefs, hoRefs;
    for (int iHcal = 0; iHcal < nHcal; ++iHcal)
      if (hcalAssignedTo[iHcal] == iEcal) hcalRefs.push_back(reco::PFClusterRef(hcalHandle, iHcal));
    for (int iHf = 0; iHf < nHf; ++iHf)
      if (hfAssignedTo[iHf] == iEcal) hfRefs.push_back(reco::PFClusterRef(hfHandle, iHf));
    for (int iHo = 0; iHo < nHo; ++iHo)
      if (hoAssignedTo[iHo] == iEcal) hoRefs.push_back(reco::PFClusterRef(hoHandle, iHo));

    reco::PFClusterRefVector ecalRefs;
    reco::PFClusterRef seedRef(ecalHandle, iEcal);
    ecalRefs.push_back(seedRef);

    double eEcal = ecalClus.energy();
    double eHcal = 0., eHf = 0., eHo = 0.;
    for (const auto& href : hcalRefs) eHcal += href->energy();
    for (const auto& href : hfRefs) eHf += href->energy();
    for (const auto& href : hoRefs) eHo += href->energy();

    if(ecalRefs.size() + hcalRefs.size() + hfRefs.size() + hoRefs.size() <= 1) continue;
    reco::PFEHCluster ehClust(seedRef, ecalRefs, hcalRefs, hfRefs, hoRefs);
    ehClust.setRawEcalEnergy(eEcal);
    ehClust.setRawHcalEnergy(eHcal);
    ehClust.setRawHfEnergy(eHf);
    ehClust.setRawHoEnergy(eHo);
    ehClust.setPosition(energyWeightedPosition(ecalRefs, hcalRefs, hfRefs, hoRefs));

    output->push_back(std::move(ehClust));
    ecalUsed[iEcal] = true;

    LogDebug("PFEHClusterProducer")
        << "SC seed ecalE=" << eEcal
        << " GeV, nHcal=" << hcalRefs.size() << " hcalE=" << eHcal
        << " GeV, nHf=" << hfRefs.size() << " hfE=" << eHf
        << " GeV, nHo=" << hoRefs.size() << " hoE=" << eHo
        << " GeV, eta=" << output->back().eta() << " phi=" << output->back().phi();
  }

  // Promote unmatched HCAL clusters.
  for (int iHcal = 0; iHcal < nHcal; ++iHcal) {
    if (hcalAssignedTo[iHcal] >= 0) continue;
    const auto& hcalClus = hcalClusters[iHcal];

    reco::PFClusterRef hRef(hcalHandle, iHcal);
    reco::PFClusterRefVector emptyEcal, hcalRefs, emptyHf, emptyHo;
    hcalRefs.push_back(hRef);

    reco::PFEHCluster ehClust(hRef, emptyEcal, hcalRefs, emptyHf, emptyHo);
    ehClust.setRawEcalEnergy(0.);
    ehClust.setRawHcalEnergy(hcalClus.energy());
    ehClust.setPosition(math::XYZPoint(hcalClus.position().x(),
                                        hcalClus.position().y(),
                                        hcalClus.position().z()));
    output->push_back(std::move(ehClust));
// Jenn: Fix this, but just for testing
    //unmergedHcalOut->push_back(hcalClusters[iHcal]);
  }

  // Promote unmatched HF clusters.
  for (int iHf = 0; iHf < nHf; ++iHf) {
    if (hfAssignedTo[iHf] >= 0) continue;
/*
    const auto& hfClus = hfClusters[iHf];

    reco::PFClusterRef hRef(hfHandle, iHf);
    reco::PFClusterRefVector emptyEcal, emptyHcal, hfRefs, emptyHo;
    hfRefs.push_back(hRef);

    reco::PFEHCluster ehClust(hRef, emptyEcal, emptyHcal, hfRefs, emptyHo);
    ehClust.setRawHfEnergy(hfClus.energy());
    ehClust.setPosition(math::XYZPoint(hfClus.position().x(),
                                        hfClus.position().y(),
                                        hfClus.position().z()));
    output->push_back(std::move(ehClust));
*/

    unmergedHfOut->push_back(hfClusters[iHf]);
  }

  // Promote unmatched HO clusters.
  for (int iHo = 0; iHo < nHo; ++iHo) {
    if (hoAssignedTo[iHo] >= 0) continue;
/*
    const auto& hoClus = hoClusters[iHo];

    reco::PFClusterRef hRef(hoHandle, iHo);
    reco::PFClusterRefVector emptyEcal, emptyHcal, emptyHf, hoRefs;
    hoRefs.push_back(hRef);

    reco::PFEHCluster ehClust(hRef, emptyEcal, emptyHcal, emptyHf, hoRefs);
    ehClust.setRawHoEnergy(hoClus.energy());
    ehClust.setPosition(math::XYZPoint(hoClus.position().x(),
                                        hoClus.position().y(),
                                        hoClus.position().z()));
    output->push_back(std::move(ehClust));
*/
    unmergedHoOut->push_back(hoClusters[iHo]);
  }

  // Keep unmatched ECAL-only clusters.
  for (int iEcal = 0; iEcal < nEcal; ++iEcal) {
    if (ecalUsed[iEcal]) continue;
/*
    const auto& ecalClus = ecalClusters[iEcal];

    reco::PFClusterRef eRef(ecalHandle, iEcal);
    reco::PFClusterRefVector ecalRefs;
    ecalRefs.push_back(eRef);
    reco::PFClusterRefVector emptyHcal, emptyHf, emptyHo;

    reco::PFEHCluster ehClust(eRef, ecalRefs, emptyHcal, emptyHf, emptyHo);
    ehClust.setRawEcalEnergy(ecalClus.energy());
    ehClust.setPosition(math::XYZPoint(ecalClus.position().x(),
                                        ecalClus.position().y(),
                                        ecalClus.position().z()));
    output->push_back(std::move(ehClust));
*/
    unmergedEcalOut->push_back(ecalClusters[iEcal]);
  }

  LogDebug("PFEHClusterProducer")
      << "Produced " << output->size() << " PFEHClusters from "
      << nEcal << " ECAL, " << nHcal << " HCAL, "
      << nHf << " HF, and " << nHo << " HO PFClusters.";

  iEvent.put(pfehClusterToken_, std::move(output));
  iEvent.put(unmergedEcalToken_, std::move(unmergedEcalOut));
  iEvent.put(unmergedHcalToken_, std::move(unmergedHcalOut));
  iEvent.put(unmergedHfToken_, std::move(unmergedHfOut));
  iEvent.put(unmergedHoToken_, std::move(unmergedHoOut));

}

// ============================================================================
// energyWeightedPosition
// ============================================================================
math::XYZPoint PFEHClusterProducer::energyWeightedPosition(
    const reco::PFClusterRefVector& ecalClusters,
    const reco::PFClusterRefVector& hcalClusters,
    const reco::PFClusterRefVector& hfClusters,
    const reco::PFClusterRefVector& hoClusters) {
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
  accumulate(hfClusters);
  accumulate(hoClusters);

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
  desc.add<edm::InputTag>("hfClusters", edm::InputTag("particleFlowClusterHF"))
      ->setComment("HF PFCluster collection");
  desc.add<edm::InputTag>("hoClusters", edm::InputTag("particleFlowClusterHO"))
      ->setComment("HO PFCluster collection");
  desc.add<double>("matchingDeltaR", 0.15)
      ->setComment("Maximum dR between ECAL and HCAL cluster positions for matching");
  desc.add<double>("matchingDeltaRHF", 0.3)
      ->setComment("Maximum dR between ECAL and HF cluster positions for matching");
  desc.add<double>("matchingDeltaRHO", 0.2)
      ->setComment("Maximum dR between ECAL and HO cluster positions for matching");
  descriptions.add("pfEHClusterProducer", desc);
}

// ============================================================================
// Plugin registration
// ============================================================================
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFEHClusterProducer);
