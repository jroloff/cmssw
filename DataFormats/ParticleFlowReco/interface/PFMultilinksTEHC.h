#ifndef __PFMultilinksTEHC__
#define __PFMultilinksTEHC__

// Done by Glowinski & Gouzevitch

#include <vector>
#include "DataFormats/ParticleFlowReco/interface/PFRecTrackFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFEHClusterFwd.h"

namespace reco {

  /// \brief Abstract This class is used by the KDTree Track / EHCluster
  /// linker to store all found links.
  ///
  struct PFEHMultilink {
    PFEHMultilink(const reco::PFEHClusterRef& clusterref) : trackRef(), clusterRef(clusterref) {}
    PFEHMultilink(const reco::PFRecTrackRef& trackref) : trackRef(trackref), clusterRef() {}
    reco::PFRecTrackRef trackRef;
    reco::PFEHClusterRef clusterRef;
  };
  /// collection of PFSuperCluster objects
  typedef std::vector<PFEHMultilink> PFEHMultilinksType;
  class PFMultiLinksTEHC {
  public:
    bool isValid;
    PFEHMultilinksType linkedPFObjects;

  public:
    PFMultiLinksTEHC(bool isvalid = false) : isValid(isvalid) {}
  };
}  // namespace reco

#endif
