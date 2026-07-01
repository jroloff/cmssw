#ifndef ParticleFlowReco_PFEHClusterFwd_h
#define ParticleFlowReco_PFEHClusterFwd_h
#include <vector>
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/RefProd.h"

namespace reco {
  class PFEHCluster;

  /// collection of PFEHCluster objects
  typedef std::vector<PFEHCluster> PFEHClusterCollection;

  /// persistent reference to PFEHCluster objects
  typedef edm::Ref<PFEHClusterCollection> PFEHClusterRef;

  /// reference to PFEHCluster collection
  typedef edm::RefProd<PFEHClusterCollection> PFEHClusterRefProd;

  /// vector of references to PFEHCluster objects all in the same collection
  typedef edm::RefVector<PFEHClusterCollection> PFEHClusterRefVector;

  /// iterator over a vector of references to PFEHCluster objects
  typedef PFEHClusterRefVector::iterator PFEHCluster_iterator;
}  // namespace reco

#endif
