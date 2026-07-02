#include "RecoParticleFlow/PFProducer/interface/BlockElementLinkerBase.h"
#include "DataFormats/ParticleFlowReco/interface/PFEHCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowReco/interface/PFTrajectoryPoint.h"
#include "RecoParticleFlow/PFClusterTools/interface/LinkByRecHit.h"

class TrackAndEHClusterLinker : public BlockElementLinkerBase {
public:
  TrackAndEHClusterLinker(const edm::ParameterSet& conf)
      : BlockElementLinkerBase(conf),
        useKDTree_(conf.getParameter<bool>("useKDTree")),
        debug_(conf.getUntrackedParameter<bool>("debug", false)) {}

  bool linkPrefilter(const reco::PFBlockElement*, const reco::PFBlockElement*) const override;

  double testLinkHCal(const reco::PFBlockElement* elem1, const reco::PFBlockElement* elem2) const;
  double testLinkECal(const reco::PFBlockElement* elem1, const reco::PFBlockElement* elem2) const;


  double testLink(const reco::PFBlockElement*, const reco::PFBlockElement*) const override;

private:
  const bool useKDTree_, debug_;
};

DEFINE_EDM_PLUGIN(BlockElementLinkerFactory, TrackAndEHClusterLinker, "TrackAndEHClusterLinker");

bool TrackAndEHClusterLinker::linkPrefilter(const reco::PFBlockElement* elem1, const reco::PFBlockElement* elem2) const {
  bool result = false;
  // Track-ECAL KDTree multilinks are stored to eh's elem
  switch (elem1->type()) {
    case reco::PFBlockElement::TRACK:
      result = (elem2->isMultilinksValidEH(elem1->type()) && !elem2->getMultilinksEH(elem1->type()).empty() &&
                elem1->isMultilinksValidEH(elem2->type()));
      break;
    case reco::PFBlockElement::EH:
      result = (elem1->isMultilinksValidEH(elem2->type()) && !elem1->getMultilinksEH(elem2->type()).empty() &&
                elem2->isMultilinksValidEH(elem1->type()));
    default:
      break;
  }
  return (useKDTree_ ? result : true);
}

// It would be better to use the official functions, but this seems a bit annoying for what I want to implement at the moment
double TrackAndEHClusterLinker::testLinkECal(const reco::PFBlockElement* elem1, const reco::PFBlockElement* elem2) const {
  constexpr reco::PFTrajectoryPoint::LayerType ECALShowerMax = reco::PFTrajectoryPoint::ECALShowerMax;
  const reco::PFBlockElementCluster* ecalelem(nullptr);
  const reco::PFBlockElementTrack* tkelem(nullptr);
  double dist(-1.0);
  if (elem1->type() < elem2->type()) {
    tkelem = static_cast<const reco::PFBlockElementTrack*>(elem1);
    ecalelem = static_cast<const reco::PFBlockElementCluster*>(elem2);
  } else {
    tkelem = static_cast<const reco::PFBlockElementTrack*>(elem2);
    ecalelem = static_cast<const reco::PFBlockElementCluster*>(elem1);
  }
  const reco::PFRecTrackRef& trackref = tkelem->trackRefPF();
  const reco::PFClusterRef& clusterref = ecalelem->clusterRef();
  const reco::PFCluster::REPPoint& ecalreppos = clusterref->positionREP();
  const reco::PFTrajectoryPoint& tkAtECAL = trackref->extrapolatedPoint(ECALShowerMax);

  // Check if the linking has been done using the KDTree algo
  // Glowinski & Gouzevitch
  if (useKDTree_ && ecalelem->isMultilinksValidEH(tkelem->type())) {  //KDTree Algo
    const reco::PFMultilinksType& multilinks = ecalelem->getMultilinksEH(tkelem->type());
    const double tracketa = tkAtECAL.positionREP().Eta();
    const double trackphi = tkAtECAL.positionREP().Phi();
    // Check if the link Track/Ecal exist
    reco::PFMultilinksType::const_iterator mlit = multilinks.begin();
    for (; mlit != multilinks.end(); ++mlit)
      if (mlit->trackRef == trackref)
        break;

    // If the link exist, we fill dist and linktest.
    if (mlit != multilinks.end()) {
      dist = LinkByRecHit::computeDist(ecalreppos.Eta(), ecalreppos.Phi(), tracketa, trackphi);
    }

  } else {  // Old algorithm
    if (tkAtECAL.isValid())
      dist = LinkByRecHit::testTrackAndClusterByRecHit(*trackref, *clusterref, false, debug_);
  }

  if (debug_) {
    if (dist > 0.) {
      std::cout << " Here a link has been established"
                << " between a track an Ecal with dist  " << dist << std::endl;
    } else
      std::cout << " No link found " << std::endl;
  }
  return dist;

}

double TrackAndEHClusterLinker::testLinkHCal(const reco::PFBlockElement* elem1, const reco::PFBlockElement* elem2) const {
  const reco::PFBlockElementCluster* hcalelem(nullptr);
  const reco::PFBlockElementTrack* tkelem(nullptr);
  double dist(-1.0);

  double trajectoryLayerEntrance_ = reco::PFTrajectoryPoint::HCALEntrance;
  double trajectoryLayerExit_ = reco::PFTrajectoryPoint::HCALExit;

  if (elem1->type() < elem2->type()) {
    tkelem = static_cast<const reco::PFBlockElementTrack*>(elem1);
    hcalelem = static_cast<const reco::PFBlockElementCluster*>(elem2);
  } else {
    tkelem = static_cast<const reco::PFBlockElementTrack*>(elem2);
    hcalelem = static_cast<const reco::PFBlockElementCluster*>(elem1);
  }
  const reco::PFRecTrackRef& trackref = tkelem->trackRefPF();
  const reco::PFClusterRef& clusterref = hcalelem->clusterRef();
  const reco::PFCluster::REPPoint& hcalreppos = clusterref->positionREP();
  const reco::PFTrajectoryPoint& tkAtHCALEnt = trackref->extrapolatedPoint(trajectoryLayerEntrance_);
  const reco::PFCluster::REPPoint& tkreppos = tkAtHCALEnt.positionREP();
  // Check exit point
  double dHEta = 0.;
  double dHPhi = 0.;
  double dRHCALEx = 0.;
  bool checkExit_ = true;
  if (checkExit_) {
    const reco::PFTrajectoryPoint& tkAtHCALEx = trackref->extrapolatedPoint(trajectoryLayerExit_);
    dHEta = (tkAtHCALEx.positionREP().Eta() - tkAtHCALEnt.positionREP().Eta());
    dHPhi = reco::deltaPhi(tkAtHCALEx.positionREP().Phi(), tkAtHCALEnt.positionREP().Phi());
    dRHCALEx = tkAtHCALEx.position().R();
  }
  // Check if the linking has been done using the KDTree algo
  // Glowinski & Gouzevitch
  if (useKDTree_ && tkelem->isMultilinksValidEH(hcalelem->type())) {  //KDTree Algo
    const reco::PFMultilinksType& multilinks = tkelem->getMultilinksEH(hcalelem->type());

    // Check if the link Track/Hcal exist
    reco::PFMultilinksType::const_iterator mlit = multilinks.begin();
    for (; mlit != multilinks.end(); ++mlit)
      if (mlit->clusterRef == clusterref)
        break;

    // If the link exist, we fill dist and linktest.

    if (mlit != multilinks.end()) {
      // when checkExit_ is false
      if (!checkExit_) {
        dist = LinkByRecHit::computeDist(hcalreppos.Eta(), hcalreppos.Phi(), tkreppos.Eta(), tkreppos.Phi());
      }
      // when checkExit_ is true
      else {
        //special case ! A looper  can exit the barrel inwards and hit the endcap
        //In this case calculate the distance based on the first crossing since
        //the looper will probably never make it to the endcap
        if (dRHCALEx < tkAtHCALEnt.position().R()) {
          dist = LinkByRecHit::computeDist(hcalreppos.Eta(), hcalreppos.Phi(), tkreppos.Eta(), tkreppos.Phi());
        } else {
          dist = LinkByRecHit::computeDist(
              hcalreppos.Eta(), hcalreppos.Phi(), tkreppos.Eta() + 0.1 * dHEta, tkreppos.Phi() + 0.1 * dHPhi);
        }
      }  // checkExit_
    }  // multilinks

  } else {  // Old algorithm
    if (tkAtHCALEnt.isValid())
      dist = LinkByRecHit::testTrackAndClusterByRecHit(*trackref, *clusterref, false, debug_);
  }
  return dist;
}


double TrackAndEHClusterLinker::testLink(const reco::PFBlockElement* elem1, const reco::PFBlockElement* elem2) const {
/*
  const reco::PFBlockElementCluster* ehclustelem(nullptr);
  const reco::PFBlockElementTrack* tkelem(nullptr);
  const reco::PFEHClusterRef& pfehcluster = ehclustelem->ehClusterRef();

  double dist = -1;
  
  const reco::PFClusterRefVector ecalClusters = pfehcluster->ecalClusters();
  for(auto eclus: ecalClusters){
    const reco::PFBlockElementCluster* blockElem = new reco::PFBlockElementCluster(eclus, reco::PFBlockElement::Type::ECAL);
    double newdist = testLinkEcal(tkelem, blockElem);
    if(dist < 0 || newdist < dist) dist = newdist;
  }

  const reco::PFClusterRefVector hcalClusters = pfehcluster->hcalClusters();
  for(auto hclus: hcalClusters){
    const reco::PFBlockElementCluster* blockElem = new reco::PFBlockElementCluster(hclus, reco::PFBlockElement::Type::HCAL);
    double newdist = testLinkHcal(tkelem, blockElem);
    if(dist < 0 || newdist < dist) dist = newdist;
  }

  return dist;
*/


  constexpr reco::PFTrajectoryPoint::LayerType ECALShowerMax = reco::PFTrajectoryPoint::ECALShowerMax;
  const reco::PFBlockElementCluster* ecalelem(nullptr);
  const reco::PFBlockElementTrack* tkelem(nullptr);
  double dist(-1.0);
  if (elem1->type() < elem2->type()) {
    tkelem = static_cast<const reco::PFBlockElementTrack*>(elem1);
    ecalelem = static_cast<const reco::PFBlockElementCluster*>(elem2);
  } else {
    tkelem = static_cast<const reco::PFBlockElementTrack*>(elem2);
    ecalelem = static_cast<const reco::PFBlockElementCluster*>(elem1);
  }
  const reco::PFRecTrackRef& trackref = tkelem->trackRefPF();
  const reco::PFClusterRef& clusterref = ecalelem->clusterRef();
  const reco::PFCluster::REPPoint& ecalreppos = clusterref->positionREP();
  const reco::PFTrajectoryPoint& tkAtECAL = trackref->extrapolatedPoint(ECALShowerMax);

  // Check if the linking has been done using the KDTree algo
  // Glowinski & Gouzevitch
  if (useKDTree_ && ecalelem->isMultilinksValidEH(tkelem->type())) {  //KDTree Algo
    const reco::PFMultilinksType& multilinks = ecalelem->getMultilinksEH(tkelem->type());
    const double tracketa = tkAtECAL.positionREP().Eta();
    const double trackphi = tkAtECAL.positionREP().Phi();
    // Check if the link Track/Ecal exist
    reco::PFMultilinksType::const_iterator mlit = multilinks.begin();
    for (; mlit != multilinks.end(); ++mlit)
      if (mlit->trackRef == trackref)
        break;

    // If the link exist, we fill dist and linktest.
    if (mlit != multilinks.end()) {
      dist = LinkByRecHit::computeDist(ecalreppos.Eta(), ecalreppos.Phi(), tracketa, trackphi);
    }

  } else {  // Old algorithm
    if (tkAtECAL.isValid())
      dist = LinkByRecHit::testTrackAndClusterByRecHit(*trackref, *clusterref, false, debug_);
  }

  if (debug_) {
    if (dist > 0.) {
      std::cout << " Here a link has been established"
                << " between a track an Ecal with dist  " << dist << std::endl;
    } else
      std::cout << " No link found " << std::endl;
  }
  return dist;

}
