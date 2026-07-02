#include "RecoParticleFlow/PFProducer/interface/BlockElementImporterBase.h"
#include "DataFormats/ParticleFlowReco/interface/PFEHCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementEHCluster.h"

class EHClusterImporter : public BlockElementImporterBase {
public:
  EHClusterImporter(const edm::ParameterSet& conf, edm::ConsumesCollector& cc)
      : BlockElementImporterBase(conf, cc),
        _src(cc.consumes<reco::PFEHClusterCollection>(conf.getParameter<edm::InputTag>("source"))) {}

  void importToBlock(const edm::Event&, ElementList&) const override;

private:
  edm::EDGetTokenT<reco::PFEHClusterCollection> _src;
};

DEFINE_EDM_PLUGIN(BlockElementImporterFactory, EHClusterImporter, "EHClusterImporter");

void EHClusterImporter::importToBlock(const edm::Event& e, BlockElementImporterBase::ElementList& elems) const {
  auto clusters = e.getHandle(_src);
  auto cbegin = clusters->cbegin();
  auto cend = clusters->cend();
  for (auto clus = cbegin; clus != cend; ++clus) {
    reco::PFBlockElement::Type type = reco::PFBlockElement::NONE;
    reco::PFEHClusterRef cref(clusters, std::distance(cbegin, clus));
    reco::PFBlockElement* cptr = new reco::PFBlockElementEHCluster(cref);
    elems.emplace_back(cptr);
  }
}
