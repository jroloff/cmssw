#include "RecoParticleFlow/PFProducer/interface/BlockElementLinkerBase.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFEHCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementEHCluster.h"
#include "RecoParticleFlow/PFClusterTools/interface/LinkByRecHit.h"

class EHAndHOLinker : public BlockElementLinkerBase {
public:
  EHAndHOLinker(const edm::ParameterSet& conf)
      : BlockElementLinkerBase(conf),
        useKDTree_(conf.getParameter<bool>("useKDTree")),
        debug_(conf.getUntrackedParameter<bool>("debug", false)) {}

  double testLink(const reco::PFBlockElement*, const reco::PFBlockElement*) const override;

private:
  bool useKDTree_, debug_;
};

DEFINE_EDM_PLUGIN(BlockElementLinkerFactory, EHAndHOLinker, "EHAndHOLinker");

double EHAndHOLinker::testLink(const reco::PFBlockElement* elem1, const reco::PFBlockElement* elem2) const {
  const reco::PFBlockElementEHCluster *hcalelem(nullptr);
  const reco::PFBlockElementCluster *hoelem(nullptr);
  double dist(-1.0);
  if (elem1->type() > elem2->type()) {
    hcalelem = static_cast<const reco::PFBlockElementEHCluster*>(elem1);
    hoelem = static_cast<const reco::PFBlockElementCluster*>(elem2);
  } else {
    hcalelem = static_cast<const reco::PFBlockElementEHCluster*>(elem2);
    hoelem = static_cast<const reco::PFBlockElementCluster*>(elem1);
  }
  const reco::PFEHClusterRef& ehref = hcalelem->ehClusterRef();
  const reco::PFClusterRef& horef = hoelem->clusterRef();

       for (auto const& hcalref : ehref->hcalClusters()) {

  const reco::PFCluster::REPPoint& hcalreppos = hcalref->positionREP();
  if (hcalref.isNull() || horef.isNull()) {
    throw cms::Exception("BadClusterRefs") << "PFBlockElementCluster's refs are null!";
  }
  dist = std::min(dist,(std::abs(hcalreppos.Eta()) < 1.5
              ? LinkByRecHit::computeDist(
                    hcalreppos.Eta(), hcalreppos.Phi(), horef->positionREP().Eta(), horef->positionREP().Phi())
              : -1.0));
        }
  return (dist < 0.2 ? dist : -1.0);
}
