#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFEHCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "RecoParticleFlow/PFProducer/interface/KDTreeLinkerBase.h"
#include "CommonTools/RecoAlgos/interface/KDTreeLinkerAlgo.h"

#include "TMath.h"

// This class is used to find all links between Tracks and EH clusters
// using a KDTree algorithm.
// It is used in PFBlockAlgo.cc in the function links().
class KDTreeLinkerTrackEHCluster : public KDTreeLinkerBase {
public:
  KDTreeLinkerTrackEHCluster(const edm::ParameterSet &conf);
  ~KDTreeLinkerTrackEHCluster() override;

  // With this method, we create the list of track that we want to link.
  void insertTargetElt(reco::PFBlockElement *track) override;

  // Here, we create the list of EHCluster that we want to link. From ecalCluster
  // and fraction, we will create a second list of rechits that will be used to
  // build the KDTree.
  void insertFieldClusterElt(reco::PFBlockElement *ehCluster) override;

  // The KDTree building from rechits list.
  void buildTree() override;

  // Here we will iterate over all tracks. For each track intersection point with EH,
  // we will search the closest rechits in the KDTree, from rechits we will find the
  // ehClusters and after that we will check the links between the track and
  // all closest ehClusters.
  void searchLinks() override;

  // Here, we will store all Track/EH founded links in the PFBlockElement class
  // of each psCluster in the PFmultilinks field.
  void updatePFBlockEltWithLinks() override;

  // Here we free all allocated structures.
  void clear() override;

private:
  // Data used by the KDTree algorithm : sets of Tracks and EH clusters.
  BlockEltSet targetSet_;
  BlockEltSet fieldClusterSet_;

  // Sets of rechits that compose the EH clusters.
  RecHitSet rechitsSetEcal_;
  RecHitSet rechitsSetHcal_;

  // TrajectoryPoints
  std::string trajectoryLayerEntranceString_;
  std::string trajectoryLayerExitString_;
  reco::PFTrajectoryPoint::LayerType trajectoryLayerEntrance_;
  reco::PFTrajectoryPoint::LayerType trajectoryLayerExit_;
  bool checkExit_;

  // Map of linked Track/EH clusters.
  BlockElt2BlockEltMap cluster2TargetLinks_;

  // Map of the EH clusters associated to a rechit.
  RecHit2BlockEltMap rechit2ClusterLinks_;

  // KD trees
  KDTreeLinkerAlgo<reco::PFRecHit const *> treeEcal_;
  KDTreeLinkerAlgo<reco::PFRecHit const *> treeHcal_;
};

// the text name is different so that we can easily
// construct it when calling the factory
DEFINE_EDM_PLUGIN(KDTreeLinkerFactory, KDTreeLinkerTrackEHCluster, "KDTreeTrackAndEHClusterLinker");

KDTreeLinkerTrackEHCluster::KDTreeLinkerTrackEHCluster(const edm::ParameterSet &conf) : KDTreeLinkerBase(conf) {}

KDTreeLinkerTrackEHCluster::~KDTreeLinkerTrackEHCluster() { clear(); }

void KDTreeLinkerTrackEHCluster::insertTargetElt(reco::PFBlockElement *track) {
  if (track->trackRefPF()->extrapolatedPoint(reco::PFTrajectoryPoint::ECALShowerMax).isValid()) {
    targetSet_.insert(track);
  }
}

void KDTreeLinkerTrackEHCluster::insertFieldClusterElt(reco::PFBlockElement *ehCluster) {
  // Make separate links for ECal and HCal since the track extrapolation is different
  const reco::PFClusterRefVector ecalClusters = ehCluster->ehClusterRef()->ecalClusters();
  // We create a list of ehCluster
  fieldClusterSet_.insert(ehCluster);
  for(auto clusterref: ecalClusters){
    const std::vector<reco::PFRecHitFraction> &fraction = clusterref->recHitFractions();

    //fieldClusterSet_.insert(ehCluster);
    for (size_t rhit = 0; rhit < fraction.size(); ++rhit) {
      const reco::PFRecHitRef &rh = fraction[rhit].recHitRef();
      double fract = fraction[rhit].fraction();
  
      // TODO: May want this to depend on the EH fraction, not the ECAL fraction
      if ((rh.isNull()) || (fract < cutOffFrac)){
        continue;
      }
  
      const reco::PFRecHit &rechit = *rh;
  
      // We save the links rechit to EcalClusters
      rechit2ClusterLinks_[&rechit].insert(ehCluster);
  
      // We create a list of rechits
      rechitsSetEcal_.insert(&rechit);
    }
  }

  const reco::PFClusterRefVector hcalClusters = ehCluster->ehClusterRef()->hcalClusters();
  for(auto clusterref: hcalClusters){
    const std::vector<reco::PFRecHitFraction> &fraction = clusterref->recHitFractions();

    // We create a list of hcalCluster
    //fieldClusterSet_.insert(ehCluster);
    for (size_t rhit = 0; rhit < fraction.size(); ++rhit) {
      const reco::PFRecHitRef &rh = fraction[rhit].recHitRef();
      double fract = fraction[rhit].fraction();
 
      if ((rh.isNull()) || (fract < cutOffFrac)){
        continue;
      }
      const reco::PFRecHit &rechit = *rh;
 
      // We save the links rechit to EcalClusters
      rechit2ClusterLinks_[&rechit].insert(ehCluster);
 
      // We create a list of rechits
      rechitsSetHcal_.insert(&rechit);
    }
  }

}

void KDTreeLinkerTrackEHCluster::buildTree() {
  // List of pseudo-rechits that will be used to create the KDTree
  std::vector<KDTreeNodeInfo<reco::PFRecHit const *, 2>> eltListEcal;

  // Here we define the upper/lower bounds of the 2D space (eta/phi).
  float phimin = -1.0 * M_PI - phiOffset_;
  float phimax = M_PI + phiOffset_;

  // etamin-etamax, phimin-phimax
  KDTreeBox region(-3.0f, 3.0f, phimin, phimax);

  // Filling of this list
  for (RecHitSet::const_iterator it = rechitsSetEcal_.begin(); it != rechitsSetEcal_.end(); it++) {
    const reco::PFRecHit::REPPoint &posrep = (*it)->positionREP();

    KDTreeNodeInfo<reco::PFRecHit const *, 2> rh1(*it, posrep.eta(), posrep.phi());
    eltListEcal.push_back(rh1);

    // Here we solve the problem of phi circular set by duplicating some rechits
    // too close to -Pi (or to Pi) and adding (substracting) to them 2 * Pi.
    if (rh1.dims[1] > (M_PI - phiOffset_)) {
      float phi = rh1.dims[1] - 2 * M_PI;
      KDTreeNodeInfo<reco::PFRecHit const *, 2> rh2(*it, float(posrep.eta()), phi);
      eltListEcal.push_back(rh2);
    }

    if (rh1.dims[1] < (M_PI * -1.0 + phiOffset_)) {
      float phi = rh1.dims[1] + 2 * M_PI;
      KDTreeNodeInfo<reco::PFRecHit const *, 2> rh3(*it, float(posrep.eta()), phi);
      eltListEcal.push_back(rh3);
    }
  }


  // We may now build the KDTree
  treeEcal_.build(eltListEcal, region);

  std::vector<KDTreeNodeInfo<reco::PFRecHit const*, 2>> eltListHcal;
  // Filling of this list
  for (RecHitSet::const_iterator it = rechitsSetHcal_.begin(); it != rechitsSetHcal_.end(); it++) {
    const reco::PFRecHit::REPPoint& posrep = (*it)->positionREP();

    KDTreeNodeInfo<reco::PFRecHit const*, 2> rh1(*it, posrep.eta(), posrep.phi());
    eltListHcal.push_back(rh1);

    // Here we solve the problem of phi circular set by duplicating some rechits
    // too close to -Pi (or to Pi) and adding (substracting) to them 2 * Pi.
    if (rh1.dims[1] > (M_PI - phiOffset_)) {
      float phi = rh1.dims[1] - 2 * M_PI;
      KDTreeNodeInfo<reco::PFRecHit const*, 2> rh2(*it, float(posrep.eta()), phi);
      eltListHcal.push_back(rh2);
    }

    if (rh1.dims[1] < (M_PI * -1.0 + phiOffset_)) {
      float phi = rh1.dims[1] + 2 * M_PI;
      KDTreeNodeInfo<reco::PFRecHit const*, 2> rh3(*it, float(posrep.eta()), phi);
      eltListHcal.push_back(rh3);
    }
  }

  // We may now build the KDTree
  treeHcal_.build(eltListHcal, region);
}

void KDTreeLinkerTrackEHCluster::searchLinks() {
  // Most of the code has been taken from LinkByRecHit.cc

  // We iterate over the tracks.
  for (BlockEltSet::iterator it = targetSet_.begin(); it != targetSet_.end(); it++) {
    reco::PFRecTrackRef trackref = (*it)->trackRefPF();

    const reco::PFTrajectoryPoint &atECAL = trackref->extrapolatedPoint(reco::PFTrajectoryPoint::ECALShowerMax);

    // The track didn't reach ecal
    if (!atECAL.isValid())
      continue;

    const reco::PFTrajectoryPoint &atVertex = trackref->extrapolatedPoint(reco::PFTrajectoryPoint::ClosestApproach);

    double trackPt = sqrt(atVertex.momentum().Vect().Perp2());
    float tracketa = atECAL.positionREP().eta();
    float trackphi = atECAL.positionREP().phi();
    double trackx = atECAL.position().X();
    double tracky = atECAL.position().Y();
    double trackz = atECAL.position().Z();

    // Estimate the maximal envelope in phi/eta that will be used to find rechit candidates.
    // Same envelope for cap et barrel rechits.
    float range = cristalPhiEtaMaxSize_ * (2.0 + 1.0 / std::min(1., trackPt / 2.));

    // We search for all candidate recHits, ie all recHits contained in the maximal size envelope.
    std::vector<reco::PFRecHit const *> recHits;
    KDTreeBox trackBox(tracketa - range, tracketa + range, trackphi - range, trackphi + range);
    treeEcal_.search(trackBox, recHits);

    // Here we check all rechit candidates using the non-approximated method.
    // TODO Need to fix for HCal
    for (auto const &recHit : recHits) {
      const auto &cornersxyz = recHit->getCornersXYZ();
      const auto &posxyz = recHit->position();
      const auto &rhrep = recHit->positionREP();
      const auto &corners = recHit->getCornersREP();

      double rhsizeeta = fabs(corners[3].eta() - corners[1].eta());
      double rhsizephi = fabs(corners[3].phi() - corners[1].phi());
      if (rhsizephi > M_PI)
        rhsizephi = 2. * M_PI - rhsizephi;

      double deta = fabs(rhrep.eta() - tracketa);
      double dphi = fabs(rhrep.phi() - trackphi);
      if (dphi > M_PI)
        dphi = 2. * M_PI - dphi;

      // Find all clusters associated to given rechit
      RecHit2BlockEltMap::iterator ret = rechit2ClusterLinks_.find(recHit);

      for (BlockEltSet::const_iterator clusterIt = ret->second.begin(); clusterIt != ret->second.end(); clusterIt++) {
        const reco::PFEHClusterRef ehclusterref = (*clusterIt)->ehClusterRef();
        const reco::PFClusterRefVector ecalClusters = ehclusterref->ecalClusters();
        for(auto clusterref: ecalClusters){
        double clusterz = clusterref->position().z();
        int fracsNbr = clusterref->recHitFractions().size();

        if (clusterref->layer() == PFLayer::ECAL_BARREL) {  // BARREL
          // Check if the track is in the barrel
          if (fabs(trackz) > 300.)
            continue;

          double _rhsizeeta = rhsizeeta * (2.00 + 1.0 / (fracsNbr * std::min(1., trackPt / 2.)));
          double _rhsizephi = rhsizephi * (2.00 + 1.0 / (fracsNbr * std::min(1., trackPt / 2.)));

          // Check if the track and the cluster are linked
          if (deta < (_rhsizeeta / 2.) && dphi < (_rhsizephi / 2.))
            cluster2TargetLinks_[*clusterIt].insert(*it);

        } else {  // ENDCAP

          // Check if the track is in the cap
          if (fabs(trackz) < 300.)
            continue;
          if (trackz * clusterz < 0.)
            continue;

          double x[5];
          double y[5];
          for (unsigned jc = 0; jc < 4; ++jc) {
            auto cornerposxyz = cornersxyz[jc];
            x[3 - jc] = cornerposxyz.x() +
                        (cornerposxyz.x() - posxyz.x()) * (1.00 + 0.50 / fracsNbr / std::min(1., trackPt / 2.));
            y[3 - jc] = cornerposxyz.y() +
                        (cornerposxyz.y() - posxyz.y()) * (1.00 + 0.50 / fracsNbr / std::min(1., trackPt / 2.));
          }

          x[4] = x[0];
          y[4] = y[0];

          bool isinside = TMath::IsInside(trackx, tracky, 5, x, y);

          // Check if the track and the cluster are linked
          if (isinside)
            cluster2TargetLinks_[*clusterIt].insert(*it);
        }



        }
      }
    }
  }



  // We iterate over the tracks.
  for (BlockEltSet::iterator it = targetSet_.begin(); it != targetSet_.end(); it++) {
    reco::PFRecTrackRef trackref = (*it)->trackRefPF();

    const reco::PFTrajectoryPoint& atHCAL = trackref->extrapolatedPoint(trajectoryLayerEntrance_);

    // The track didn't reach hcal
    if (!atHCAL.isValid())
      continue;

    // In case the exit point check is requested, check eta and phi differences between entrance and exit
    double dHeta = 0.0;
    float dHphi = 0.0;
    if (checkExit_) {
      const reco::PFTrajectoryPoint& atHCALExit = trackref->extrapolatedPoint(trajectoryLayerExit_);
      dHeta = atHCALExit.positionREP().eta() - atHCAL.positionREP().eta();
      dHphi = atHCALExit.positionREP().phi() - atHCAL.positionREP().phi();
      if (dHphi > M_PI)
        dHphi = dHphi - 2. * M_PI;
      else if (dHphi < -M_PI)
        dHphi = dHphi + 2. * M_PI;
    }  // checkExit_

    float tracketa = atHCAL.positionREP().eta() + 0.1 * dHeta;
    float trackphi = atHCAL.positionREP().phi() + 0.1 * dHphi;

    if (trackphi > M_PI)
      trackphi -= 2 * M_PI;
    else if (trackphi < -M_PI)
      trackphi += 2 * M_PI;

    // Estimate the maximal envelope in phi/eta that will be used to find rechit candidates.
    // Same envelope for cap et barrel rechits.
    double inflation = 1.;
    float rangeeta = (cristalPhiEtaMaxSize_ * (1.5 + 0.5) + 0.2 * fabs(dHeta)) * inflation;
    float rangephi = (cristalPhiEtaMaxSize_ * (1.5 + 0.5) + 0.2 * fabs(dHphi)) * inflation;

    // We search for all candidate recHits, ie all recHits contained in the maximal size envelope.
    std::vector<reco::PFRecHit const*> recHits;
    KDTreeBox trackBox(tracketa - rangeeta, tracketa + rangeeta, trackphi - rangephi, trackphi + rangephi);
    treeHcal_.search(trackBox, recHits);

    // Here we check all rechit candidates using the non-approximated method.
    for (auto const& recHit : recHits) {
      const auto& rhrep = recHit->positionREP();
      const auto& corners = recHit->getCornersREP();

      double rhsizeeta = fabs(corners[3].eta() - corners[1].eta());
      double rhsizephi = fabs(corners[3].phi() - corners[1].phi());
      if (rhsizephi > M_PI)
        rhsizephi = 2. * M_PI - rhsizephi;

      double deta = fabs(rhrep.eta() - tracketa);
      double dphi = fabs(rhrep.phi() - trackphi);
      if (dphi > M_PI)
        dphi = 2. * M_PI - dphi;

      // Find all clusters associated to given rechit
      RecHit2BlockEltMap::iterator ret = rechit2ClusterLinks_.find(recHit);

      for (BlockEltSet::iterator clusterIt = ret->second.begin(); clusterIt != ret->second.end(); clusterIt++) {
        const reco::PFEHClusterRef ehclusterref = (*clusterIt)->ehClusterRef();
        const reco::PFClusterRefVector hcalClusters = ehclusterref->hcalClusters();
        for(auto clusterref: hcalClusters){
          int fracsNbr = clusterref->recHitFractions().size();
  
          double _rhsizeeta = rhsizeeta * (1.5 + 0.5 / fracsNbr) + 0.2 * fabs(dHeta);
          double _rhsizephi = rhsizephi * (1.5 + 0.5 / fracsNbr) + 0.2 * fabs(dHphi);
  
          // Check if the track and the cluster are linked
          if (deta < (_rhsizeeta / 2.) && dphi < (_rhsizephi / 2.))
            cluster2TargetLinks_[*it].insert(*clusterIt);
        }
      }
    }
  }

}

void KDTreeLinkerTrackEHCluster::updatePFBlockEltWithLinks() {
  //TODO YG : Check if cluster positionREP() is valid ?
  std::cout << __LINE__ << std::endl;

  // Here we save in each ECAL cluster the list of phi/eta values of linked tracks.
  for (BlockElt2BlockEltMap::iterator it = cluster2TargetLinks_.begin(); it != cluster2TargetLinks_.end(); ++it) {
    const auto &ecalElt = it->first;
    const auto &trackEltSet = it->second;
    reco::PFMultiLinksTEHC multitracks(true);

    for (const auto &trackElt : trackEltSet) {
      const reco::PFRecTrackRef &trackref = trackElt->trackRefPF();

      reco::PFMultilink multiLinkEH(trackref);
      multitracks.linkedPFObjects.push_back(multiLinkEH);

      // We set the multilinks flag of the track (for links to EH) to true. It will allow us to
      // use it in an optimized way in prefilter
      trackElt->setIsValidMultilinksEH(true, _fieldType);
    }

    // We set multilinks of the EH element (for links to tracks)
    ecalElt->setMultilinksEH(multitracks, _targetType);
  }
}

void KDTreeLinkerTrackEHCluster::clear() {
  std::cout << __LINE__ << std::endl;
  targetSet_.clear();
  fieldClusterSet_.clear();

  rechitsSetEcal_.clear();
  rechitsSetHcal_.clear();

  rechit2ClusterLinks_.clear();
  cluster2TargetLinks_.clear();

  treeEcal_.clear();
  treeHcal_.clear();
}
