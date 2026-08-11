// ECALAndEHClusterLinker
// --------------------------------
// Replacement for ECALAndHCALLinker now that HCAL-bearing block elements are
// PFBlockElement::EH (carrying a PFEHClusterRef) instead of
// PFBlockElement::HCAL (carrying a plain PFClusterRef). ECALAndHCALLinker is
// still registered in the linker plugin factory, but it is dead code as of
// this restructuring: it is keyed to PFBlockElement::HCAL, and no block
// element of that type is produced any more (see EHClusterImporter /
// PFEHClusterProducer), so its testLink() is simply never invoked.
//
// This plugin restores exactly the link that ECALAndHCALLinker used to
// provide -- a standalone ECAL PFBlockElement (one NOT already absorbed into
// a merged PFEHCluster; see PFEHClusterProducer) linked to a nearby
// HCAL-bearing PFEHCluster -- so that createCandidatesHCAL()'s ECAL
// "satellite" search (block.associatedElements(iHcal, ..., ECAL, ...)) finds
// it again.
//
// A PFEHCluster can carry more than one HCAL constituent (a genuine
// ECAL+multi-HCAL merge), so, mirroring the pattern used by
// TrackAndEHClusterLinker::testLink(), we test the ECAL cluster against each
// HCAL constituent individually and keep the smallest (best) distance. For
// the common "pure HCAL" PFEHCluster case (a single HCAL constituent, no
// ECAL constituent) this reduces to exactly the old ECALAndHCALLinker
// behavior.

#include "RecoParticleFlow/PFProducer/interface/BlockElementLinkerBase.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFEHCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementEHCluster.h"
#include "RecoParticleFlow/PFClusterTools/interface/LinkByRecHit.h"
#include "FWCore/Utilities/interface/Exception.h"

class ECALAndEHClusterLinker : public BlockElementLinkerBase {
public:
  ECALAndEHClusterLinker(const edm::ParameterSet& conf)
      : BlockElementLinkerBase(conf),
        minAbsEtaEcal_(conf.getParameter<double>("minAbsEtaEcal")),
        useKDTree_(conf.getParameter<bool>("useKDTree")),  // unused, kept for config compatibility with ECALAndHCALLinker
        debug_(conf.getUntrackedParameter<bool>("debug", false)) {}

  double testLink(const reco::PFBlockElement*, const reco::PFBlockElement*) const override;

private:
  double minAbsEtaEcal_;
  bool useKDTree_, debug_;
};

DEFINE_EDM_PLUGIN(BlockElementLinkerFactory, ECALAndEHClusterLinker, "ECALAndEHClusterLinker");

double ECALAndEHClusterLinker::testLink(const reco::PFBlockElement* elem1, const reco::PFBlockElement* elem2) const {
  const reco::PFBlockElementCluster* ecalelem(nullptr);
  const reco::PFBlockElementEHCluster* ehelem(nullptr);

  if (elem1->type() == reco::PFBlockElement::EH) {
    ehelem = static_cast<const reco::PFBlockElementEHCluster*>(elem1);
    ecalelem = static_cast<const reco::PFBlockElementCluster*>(elem2);
  } else {
    ehelem = static_cast<const reco::PFBlockElementEHCluster*>(elem2);
    ecalelem = static_cast<const reco::PFBlockElementCluster*>(elem1);
  }

  const reco::PFClusterRef& ecalref = ecalelem->clusterRef();
  const reco::PFEHClusterRef& ehref = ehelem->ehClusterRef();
  if (ecalref.isNull() || ehref.isNull()) {
    throw cms::Exception("BadClusterRefs") << "ECALAndEHClusterLinker: null ECAL or EH cluster ref!";
  }

  const reco::PFCluster::REPPoint& ecalreppos = ecalref->positionREP();
  //if (std::abs(ecalreppos.Eta()) <= minAbsEtaEcal_)
  if (std::abs(ecalreppos.Eta()) > minAbsEtaEcal_)
    return -1.0;

  // A pure-HCAL PFEHCluster has exactly one HCAL constituent; a genuine
  // ECAL+HCAL merge (matchingDeltaR > 0 in PFEHClusterProducer) may have
  // several. Test against each and keep the closest, same convention as
  // TrackAndEHClusterLinker::testLink().
  double dist = -1.0;
  for (const auto& hcalref : ehref->hcalClusters()) {
    if (hcalref.isNull())
      continue;
    const double d = LinkByRecHit::computeDist(
        ecalreppos.Eta(), ecalreppos.Phi(), hcalref->positionREP().Eta(), hcalref->positionREP().Phi());
    if (d >= 0. && d < 0.2 && (dist < 0. || d < dist))
      dist = d;
  }

  return dist;
}
