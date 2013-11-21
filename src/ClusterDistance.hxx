#ifndef ClusterDistance_hxx_seen
#define ClusterDistance_hxx_seen

#include <TReconCluster.hxx>

namespace CP {
    /// This finds the minimum distance between hits in the two clusters.
    /// loops through the hits in the clusters and looks to see how close the
    /// clusters actually get.  The cluster also provides the "elliptical"
    /// moments of the cluster that can be used to estimate how close the
    /// clusters get.  That method is faster, but less accurate.
    double ClusterDistance(const CP::TReconCluster& a, 
                           const CP::TReconCluster& b);
};

#endif
