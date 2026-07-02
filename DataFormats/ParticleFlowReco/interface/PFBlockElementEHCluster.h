#ifndef __PFBlockElementEHCluster__
#define __PFBlockElementEHCluster__

#include <iostream>

#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowReco/interface/PFEHClusterFwd.h"

namespace reco {

  /// \brief Cluster Element.
  ///
  /// this class contains a reference to a PFEHCluster
  class PFBlockElementEHCluster final : public PFBlockElement {
  public:
    PFBlockElementEHCluster() {}

    /// \brief constructor.
    /// \todo add a protection against the other types...
    PFBlockElementEHCluster(const PFEHClusterRef& ref)
        : PFBlockElement(reco::PFBlockElement::EH), ehClusterRef_(ref) {}

    PFBlockElement* clone() const override { return new PFBlockElementEHCluster(*this); }

    /// \return reference to the corresponding cluster
    const PFEHClusterRef& ehClusterRef() const override { return ehClusterRef_; }


    void Dump(std::ostream& out = std::cout, const char* tab = " ") const override;

  private:
    /// reference to the corresponding cluster
    PFEHClusterRef ehClusterRef_;
  };
}  // namespace reco

#endif
