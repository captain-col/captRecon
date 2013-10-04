#include "THitProximityCluster.hxx"

#include <cmath>

// Check to see if two hits are close together.  Since the hits might have
// very different extents in the XY plane and along the Z axis, this
// subtracts off the hit sizes.
double 
CP::HitProximity::Metric::operator() (const CP::HitProximity::Arg& lhs,
                                      const CP::HitProximity::Arg& rhs) {
    double xSep = 
        std::abs(lhs->GetPosition().X()-rhs->GetPosition().X());
    xSep = 
        std::max(0.0, xSep - lhs->GetRMS().X() - rhs->GetRMS().X());
    double ySep =
        std::abs(lhs->GetPosition().Y()-rhs->GetPosition().Y());
    ySep =
        std::max(0.0, ySep - lhs->GetRMS().Y() - rhs->GetRMS().Y());
    double zSep = 
        std::abs(lhs->GetPosition().Z()-rhs->GetPosition().Z());
    zSep = 
        std::max(0.0, zSep - lhs->GetRMS().Z() - rhs->GetRMS().Z());
    return std::sqrt(xSep*xSep + ySep*ySep + zSep*zSep);
}
