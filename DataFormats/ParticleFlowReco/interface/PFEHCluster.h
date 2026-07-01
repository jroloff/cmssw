#ifndef DataFormats_ParticleFlowReco_PFEHCluster_h
#define DataFormats_ParticleFlowReco_PFEHCluster_h

// PFEHCluster
// ----------------------
// A merged cluster formed by associating one or more PFClusters from the ECAL
// with one or more PFClusters from the HCAL based on proximity in eta-phi space.
//
// Design notes:
//   - Position is the energy-weighted centroid of all constituent clusters.
//   - rawEcalEnergy / rawHcalEnergy are the sums of the raw (pre-PFAlgo) cluster
//     energies from each subdetector.  The PF hadronic calibration is intentionally
//     NOT applied here; that remains the responsibility of PFAlgo when it consumes
//     these objects, so that the existing calibration chain is preserved.
//   - hOverE = rawHcalEnergy / rawEcalEnergy  (stored for convenience).
//   - The constituent cluster refs allow downstream code to retrieve shower-shape
//     variables, depth information, etc. from the original clusters.
//   - The algorithm only allows one ECal cluster, but the code is configured to
//     allow more, in case future changes to the cluster association are needed

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"

#include <vector>

namespace reco {

  class PFEHCluster : public CaloCluster {
  public:
    // ------------------------------------------------------------------ types
    using REPPoint = reco::PFCluster::REPPoint;  // (rho, eta, phi)

    // --------------------------------------------------------------- ctors
    PFEHCluster() {};

    // Build directly from constituent cluster refs.
    // The seed is the highest-energy ECAL cluster (or HCAL-only if no ECAL).
    PFEHCluster(const reco::PFClusterRef& seed,
                           const reco::PFClusterRefVector& ecalClusters,
                           const reco::PFClusterRefVector& hcalClusters){
           ecalClusters_ = ecalClusters;
           hcalClusters_ = hcalClusters;
           seed_ = seed;
    };

    // ---------------------------------------------------------------- getters

    // Energy-weighted position (eta, phi) of all constituents.
    double eta() const { return position_.eta(); }
    double phi() const { return position_.phi(); }

    // Sum of calibrated energies of all ECAL clusters.
    double rawEcalEnergy() const { return rawEcalEnergy_; }

    // Sum of calibrated energies of all HCAL clusters.
    double rawHcalEnergy() const { return rawHcalEnergy_; }

    // Total calibrated energy (ECAL + HCAL).
    double energy() const { return rawEcalEnergy_ + rawHcalEnergy_; }

    // H/E ratio (0 for ECAL-only clusters, infinity guard included).
    double hOverE() const {
      return (rawEcalEnergy_ > 0.) ? rawHcalEnergy_ / rawEcalEnergy_ : -1.;
    }

    // Number of constituent clusters in each subdetector.
    int nEcalClusters() const { return static_cast<int>(ecalClusters_.size()); }
    int nHcalClusters() const { return static_cast<int>(hcalClusters_.size()); }

    // Access to constituent cluster refs.
    const reco::PFClusterRefVector& ecalClusters() const { return ecalClusters_; }
    const reco::PFClusterRefVector& hcalClusters() const { return hcalClusters_; }

    // The seed cluster (highest-energy ECAL cluster, or highest-energy HCAL
    // cluster when no ECAL constituent is present).
    const reco::PFClusterRef& seed() const { return seed_; }

    // True when at least one cluster from each subdetector is present.
    bool isMatched() const { return !ecalClusters_.empty() && !hcalClusters_.empty(); }

    // ---------------------------------------------------------------- setters
    // (used by the producer after construction)
    void setRawEcalEnergy(double e) { rawEcalEnergy_ = e; }
    void setRawHcalEnergy(double e) { rawHcalEnergy_ = e; }
    void setPosition(const math::XYZPoint& pos) { position_ = pos; }

  private:
    math::XYZPoint          position_;       // energy-weighted centroid
    double                  rawEcalEnergy_{0.};
    double                  rawHcalEnergy_{0.};
    reco::PFClusterRef      seed_;
    reco::PFClusterRefVector ecalClusters_;
    reco::PFClusterRefVector hcalClusters_;
  };

  // Convenience typedefs expected by the framework.
  using PFEHClusterCollection = std::vector<PFEHCluster>;

}  // namespace reco

#endif  // DataFormats_ParticleFlowReco_PFEHCluster_h


