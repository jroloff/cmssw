
#include "DataFormats/ParticleFlowReco/interface/PFEHCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementEHCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowReco/interface/PFTrajectoryPoint.h"
#include "RecoParticleFlow/PFClusterTools/interface/LinkByRecHit.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "RecoParticleFlow/PFProducer/interface/BlockElementLinkerBase.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class TrackAndEHClusterLinker : public BlockElementLinkerBase {
public:
  TrackAndEHClusterLinker(const edm::ParameterSet& conf)
      : BlockElementLinkerBase(conf),
        useKDTree_(conf.getParameter<bool>("useKDTree")),
        trajectoryLayerEntranceString_(conf.getParameter<std::string>("trajectoryLayerEntrance")),
        trajectoryLayerExitString_(conf.getParameter<std::string>("trajectoryLayerExit")),
        debug_(conf.getUntrackedParameter<bool>("debug", false)) {
    // convert TrajectoryLayers info from string to enum
    trajectoryLayerEntrance_ = reco::PFTrajectoryPoint::layerTypeByName(trajectoryLayerEntranceString_);
    trajectoryLayerExit_ = reco::PFTrajectoryPoint::layerTypeByName(trajectoryLayerExitString_);
    // make sure the requested setting is supported
    assert((trajectoryLayerEntrance_ == reco::PFTrajectoryPoint::HCALEntrance &&
            trajectoryLayerExit_ == reco::PFTrajectoryPoint::HCALExit) ||
           (trajectoryLayerEntrance_ == reco::PFTrajectoryPoint::HCALEntrance &&
            trajectoryLayerExit_ == reco::PFTrajectoryPoint::Unknown) ||
           (trajectoryLayerEntrance_ == reco::PFTrajectoryPoint::VFcalEntrance &&
            trajectoryLayerExit_ == reco::PFTrajectoryPoint::Unknown));
    // flag if exit layer should be checked or not
    checkExit_ = trajectoryLayerExit_ != reco::PFTrajectoryPoint::Unknown;
  }


  bool linkPrefilter(const reco::PFBlockElement*, const reco::PFBlockElement*) const;

  // ehElem is always the parent PFBlockElementEHCluster (type EH); it carries
  // the multilinks recorded by KDTreeLinkerTrackEHCluster. clusterElem is the
  // synthetic per-constituent element (type ECAL/HCAL/HO) built from one ref
  // inside the EH cluster.
  double testLinkHCal(const reco::PFBlockElement* tkelem,
                       const reco::PFBlockElement* clusterElem,
                       const reco::PFBlockElement* ehElem) const;

  double testLink(const reco::PFBlockElement*, const reco::PFBlockElement*) const;

private:
  bool useKDTree_;
  std::string trajectoryLayerEntranceString_;
  std::string trajectoryLayerExitString_;
  reco::PFTrajectoryPoint::LayerType trajectoryLayerEntrance_;
  reco::PFTrajectoryPoint::LayerType trajectoryLayerExit_;
  bool debug_;
  bool checkExit_;
};

DEFINE_EDM_PLUGIN(BlockElementLinkerFactory, TrackAndEHClusterLinker, "TrackAndEHClusterLinker");

bool TrackAndEHClusterLinker::linkPrefilter(const reco::PFBlockElement* elem1,
                                             const reco::PFBlockElement* elem2) const {
  if (!useKDTree_) return true;

  const reco::PFBlockElement* ehElem = (elem1->type() == reco::PFBlockElement::EH) ? elem1 : elem2;
  const reco::PFBlockElement* tkElem = (elem1->type() == reco::PFBlockElement::EH) ? elem2 : elem1;

  // Multilinks are recorded EH-side (list of tracks, keyed by TRACK) and
  // track-side (validity flag only, keyed by EH) -- see
  // KDTreeLinkerTrackEHCluster::updatePFBlockEltWithLinks().
  return tkElem->isMultilinksValidEH(reco::PFBlockElement::EH) &&
         !tkElem->getMultilinksEH(reco::PFBlockElement::EH).empty() &&
         ehElem->isMultilinksValidEH(reco::PFBlockElement::TRACK);
}

double TrackAndEHClusterLinker::testLinkHCal(const reco::PFBlockElement* elem1,
                                              const reco::PFBlockElement* elem2,
                                              const reco::PFBlockElement* ehElem) const {
  const reco::PFBlockElementCluster* hcalelem(nullptr);
  const reco::PFBlockElementTrack* tkelem(nullptr);
  double dist(-1.0);

  const reco::PFTrajectoryPoint::LayerType trajectoryLayerEntrance_ = reco::PFTrajectoryPoint::HCALEntrance;
  const reco::PFTrajectoryPoint::LayerType trajectoryLayerExit_ = reco::PFTrajectoryPoint::HCALExit;

  if (elem1->type() < elem2->type()) {
    tkelem = static_cast<const reco::PFBlockElementTrack*>(elem1);
    hcalelem = static_cast<const reco::PFBlockElementCluster*>(elem2);
  } else {
    tkelem = static_cast<const reco::PFBlockElementTrack*>(elem2);
    hcalelem = static_cast<const reco::PFBlockElementCluster*>(elem1);
  }
  const reco::PFRecTrackRef& trackref = tkelem->trackRefPF();
  const reco::PFClusterRef& clusterref = hcalelem->clusterRef();
  if (trackref.isNull() || clusterref.isNull()) {
    edm::LogWarning("TrackAndEHClusterLinker") << "Null track or HCAL cluster ref; skipping.";
    return -1.;
  }
  const reco::PFCluster::REPPoint& hcalreppos = clusterref->positionREP();
  const reco::PFTrajectoryPoint& tkAtHCALEnt = trackref->extrapolatedPoint(trajectoryLayerEntrance_);
  if (!tkAtHCALEnt.isValid()) return -1.;
  const reco::PFCluster::REPPoint& tkreppos = tkAtHCALEnt.positionREP();
  // Check exit point
  double dHEta = 0.;
  double dHPhi = 0.;
  double dRHCALEx = 0.;
  const bool checkExit_ = true;
  if (checkExit_) {
    const reco::PFTrajectoryPoint& tkAtHCALEx = trackref->extrapolatedPoint(trajectoryLayerExit_);
    dHEta = (tkAtHCALEx.positionREP().Eta() - tkAtHCALEnt.positionREP().Eta());
    dHPhi = reco::deltaPhi(tkAtHCALEx.positionREP().Phi(), tkAtHCALEnt.positionREP().Phi());
    dRHCALEx = tkAtHCALEx.position().R();
  }
  // Check if the linking has been done using the KDTree algo
  // Glowinski & Gouzevitch

 if (useKDTree_ && tkelem->isMultilinksValidEH(reco::PFBlockElement::EH)) {  //KDTree Algo
    const reco::PFEHClusterRef& ehref = static_cast<const reco::PFBlockElementEHCluster*>(ehElem)->ehClusterRef();
    const reco::PFEHMultilinksType& multilinks = tkelem->getMultilinksEH(reco::PFBlockElement::EH);

    // Check if the link Track/Hcal exist
    reco::PFEHMultilinksType::const_iterator mlit = multilinks.begin();
    for (; mlit != multilinks.end(); ++mlit)
      //if (mlit->trackRef == trackref)
      if (mlit->clusterRef == ehref)
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
    dist = LinkByRecHit::testTrackAndClusterByRecHit(*trackref, *clusterref, false, debug_);
  }
  return dist;
}

double TrackAndEHClusterLinker::testLink(const reco::PFBlockElement* elem1, const reco::PFBlockElement* elem2) const {
  const reco::PFBlockElementEHCluster* ehclustelem(nullptr);
  const reco::PFBlockElementTrack* tkelem(nullptr);

  if (elem1->type() < elem2->type()) {
    tkelem = static_cast<const reco::PFBlockElementTrack*>(elem1);
    ehclustelem = static_cast<const reco::PFBlockElementEHCluster*>(elem2);
  } else {
    tkelem = static_cast<const reco::PFBlockElementTrack*>(elem2);
    ehclustelem = static_cast<const reco::PFBlockElementEHCluster*>(elem1);
  }

  const reco::PFEHClusterRef& pfehcluster = ehclustelem->ehClusterRef();
  if (pfehcluster.isNull() || tkelem->trackRefPF().isNull()) {
    edm::LogWarning("TrackAndEHClusterLinker")
        << "Null PFEHClusterRef or track ref on block element; skipping this pair.";
    return -1.;
  }

  double dist = -1.;

  for (auto hclus : pfehcluster->hcalClusters()) {
    if (hclus.isNull()) continue;
    reco::PFBlockElementCluster hcalBlockElem(hclus, reco::PFBlockElement::HCAL);
    double newdist = testLinkHCal(tkelem, &hcalBlockElem, ehclustelem);
    if (newdist >= 0. && (dist < 0. || newdist < dist)) dist = newdist;
  }

  return dist;
}
