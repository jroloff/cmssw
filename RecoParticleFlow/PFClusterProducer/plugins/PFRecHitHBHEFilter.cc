#include <memory>
#include <vector>
#include <cmath>
#include <algorithm>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFLayer.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/DetId/interface/DetId.h"

class PFRecHitHBHEFilter : public edm::stream::EDProducer<> {
public:
  explicit PFRecHitHBHEFilter(const edm::ParameterSet&);
  ~PFRecHitHBHEFilter() override = default;

  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  bool passLayer(PFLayer::Layer layer) const;
  bool passDetId(unsigned rawDetId) const;
  bool passTime(const reco::PFRecHit& hit) const;

  edm::EDGetTokenT<reco::PFRecHitCollection> srcToken_;

  double minEnergy_;
  bool useAbsEnergy_;

  bool filterLayer_;
  int requiredLayer_;

  bool filterTime_;
  double minTime_;
  double maxTime_;

  bool filterIeta_;
  int minAbsIeta_;
  int maxAbsIeta_;

  bool filterIphi_;
  std::vector<int> allowedIphi_;

  bool filterDepth_;
  std::vector<int> allowedDepth_;

  bool verbose_;
};

PFRecHitHBHEFilter::PFRecHitHBHEFilter(const edm::ParameterSet& cfg)
    : srcToken_(consumes<reco::PFRecHitCollection>(cfg.getParameter<edm::InputTag>("src"))),
      minEnergy_(cfg.getParameter<double>("minEnergy")),
      useAbsEnergy_(cfg.getParameter<bool>("useAbsEnergy")),
      filterLayer_(cfg.getParameter<bool>("filterLayer")),
      requiredLayer_(cfg.getParameter<int>("requiredLayer")),
      filterTime_(cfg.getParameter<bool>("filterTime")),
      minTime_(cfg.getParameter<double>("minTime")),
      maxTime_(cfg.getParameter<double>("maxTime")),
      filterIeta_(cfg.getParameter<bool>("filterIeta")),
      minAbsIeta_(cfg.getParameter<int>("minAbsIeta")),
      maxAbsIeta_(cfg.getParameter<int>("maxAbsIeta")),
      filterIphi_(cfg.getParameter<bool>("filterIphi")),
      allowedIphi_(cfg.getParameter<std::vector<int>>("allowedIphi")),
      filterDepth_(cfg.getParameter<bool>("filterDepth")),
      allowedDepth_(cfg.getParameter<std::vector<int>>("allowedDepth")),
      verbose_(cfg.getParameter<bool>("verbose")) {
  produces<reco::PFRecHitCollection>();
}

bool PFRecHitHBHEFilter::passLayer(PFLayer::Layer layer) const {
  if (!filterLayer_) return true;
  return static_cast<int>(layer) == requiredLayer_;
}

bool PFRecHitHBHEFilter::passTime(const reco::PFRecHit& hit) const {
  if (!filterTime_) return true;
  const double t = hit.time();
  return (t >= minTime_ && t <= maxTime_);
}

bool PFRecHitHBHEFilter::passDetId(unsigned rawDetId) const {
  DetId did(rawDetId);
  if (did.det() != DetId::Hcal) return false;

  HcalDetId hid(rawDetId);

  if (filterIeta_) {
    const int absIeta = std::abs(hid.ieta());
    if (absIeta < minAbsIeta_ || absIeta > maxAbsIeta_) return false;
  }

  if (filterIphi_) {
    if (std::find(allowedIphi_.begin(), allowedIphi_.end(), hid.iphi()) == allowedIphi_.end()) {
      return false;
    }
  }

  if (filterDepth_) {
    if (std::find(allowedDepth_.begin(), allowedDepth_.end(), hid.depth()) == allowedDepth_.end()) {
      return false;
    }
  }

  return true;
}

void PFRecHitHBHEFilter::produce(edm::Event& evt, const edm::EventSetup&) {

  edm::Handle<reco::PFRecHitCollection> handle;
  evt.getByToken(srcToken_, handle);
  if (!handle.isValid()) {
    evt.put(std::make_unique<reco::PFRecHitCollection>());
    return;
  }
  auto const& in = *handle;

  auto out = std::make_unique<reco::PFRecHitCollection>();
  std::vector<int> indexMap(in.size(), -1);

  //for (auto const& hit : in) {
  for (size_t i = 0; i < in.size(); i++) {
    auto const& hit = in[i];
    const double e = useAbsEnergy_ ? std::abs(hit.energy()) : hit.energy();
    if (e < minEnergy_) continue;

    if (!passLayer(hit.layer())) continue;
    if (!passTime(hit)) continue;
    if (!passDetId(hit.detId())) continue;

    reco::PFRecHit hit_copy = reco::PFRecHit(hit);
    hit_copy.clearNeighbours();
    hit_copy.neighbours(); // This doesn't actually seem to work for the moment. I'm not sure if this matters yet for Run 3...
    indexMap[i] = out->size();
    out->push_back(hit_copy);
  }

  // Build detId -> new index lookup
  std::unordered_map<unsigned int, unsigned int> detIdToIndex;
  for (size_t i = 0; i < out->size(); i++) {
    detIdToIndex[(*out)[i].detId()] = i;
  }


for (size_t i = 0; i < in.size(); i++) {
    if (indexMap[i] < 0) continue;
    reco::PFRecHit& oldHit = const_cast<reco::PFRecHit&>(in[i]);  // needed for non-const neighbourInfos()
    auto& newHit = (*out)[indexMap[i]];

    auto const& nbrs = oldHit.neighbours();
    auto const& infos = oldHit.neighbourInfos();

    size_t n = 0;
    for (auto const& nbrIdx : nbrs) {
       int detId = in[nbrIdx].detId();
        if (nbrIdx < in.size()) {
            auto it = detIdToIndex.find(detId);
            if (it != detIdToIndex.end()) {
                unsigned short flag = infos[n];
                short absx = (flag >> 1) & 0x7;
                short x = (flag & 1) ? absx : -absx;
                short absy = (flag >> 5) & 0x7;
                short y = (flag & (1 << 4)) ? absy : -absy;
                short absz = (flag >> 9) & 0x7;
                short z = (flag & (1 << 8)) ? absz : -absz;
                newHit.addNeighbour(x, y, z, it->second);
            }
        }
        n++;
    }
}


  evt.put(std::move(out));

}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFRecHitHBHEFilter);





