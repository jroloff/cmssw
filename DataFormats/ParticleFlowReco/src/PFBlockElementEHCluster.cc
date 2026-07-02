#include "DataFormats/ParticleFlowReco/interface/PFBlockElementEHCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFEHCluster.h"
#include "DataFormats/Common/interface/Ref.h"
#include "Math/Vector3D.h"

#include <iomanip>

using namespace reco;
using namespace std;

void PFBlockElementEHCluster::Dump(ostream& out, const char* tab) const {
  if (!out)
    return;
  // need to convert the math::XYZPoint data member of the PFEHCluster class=
  // to a displacement vector:
  ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> clusterPos(
      ehClusterRef_->position().X(), ehClusterRef_->position().Y(), ehClusterRef_->position().Z());

  clusterPos = clusterPos.Unit();
  double E = ehClusterRef_->energy();
  clusterPos *= E;
  double ET = sqrt(clusterPos.X() * clusterPos.X() + clusterPos.Y() * clusterPos.Y());

  out << setprecision(3);
  out << setiosflags(ios::right);
  out << setiosflags(ios::fixed);
  out << setw(4) << ", ET =" << setw(7) << ET;
  out << setw(4) << ", E =" << setw(7) << E;
  out << " (eta,phi,z)= (";
  out << ehClusterRef_->position().Eta() << ",";
  out << ehClusterRef_->position().Phi() << ",";
  out << ehClusterRef_->position().Z() << ")";
  out << resetiosflags(ios::right | ios::fixed);
}
