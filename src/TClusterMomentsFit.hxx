#ifndef TClusterMomentsFit_hxx
#define TClusterMomentsFit_hxx

#include <TReconCluster.hxx>

#include <TVector3.h>
#include <TLorentzVector.h>
#include <utility>

namespace CP {
    class TClusterMomentsFit;
}

/// Take a cluster, and determine the "best" parameters for the 3D line fit to
/// the cluster.  This only make sense when the cluster is "long and skinny",
/// but has the advantage that it is completely independent of the axis of the
/// cluster.  This method doesn't care about the orientation of the line.
///
/// \bug An estimate of the covariance should be made for the results of
/// TClusterMomentsFit.
class CP::TClusterMomentsFit {
public:
    /// Constructed with the cluster to be fit.
    explicit TClusterMomentsFit(const CP::TReconCluster& cluster);

    /// Get the position of the best fit line.  The position is at the center
    /// of the cluster.
    const TLorentzVector& GetPosition() const {return fPosition;}

    /// Get the direction of the best fit line.  The direction is
    /// preferentially in the positive X direction (then Y, then Z)
    const TVector3& GetDirection() const {return fDirection;}
    
    /// Get the ends of the best fit line.  The end is defined as the maximum
    /// extent of the hits along the fit direction.
    std::pair<TLorentzVector,TLorentzVector>& GetEndPoints() {
        return fEnds;
    }

private:
    /// Cache of the fitted position.
    TLorentzVector fPosition;
    /// Cache of the fitted direction.
    TVector3 fDirection;
    /// Cache of the end points.
    std::pair<TLorentzVector, TLorentzVector> fEnds;
};
#endif
