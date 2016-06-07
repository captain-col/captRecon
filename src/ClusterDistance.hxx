#ifndef ClusterDistance_hxx_seen
#define ClusterDistance_hxx_seen

#include <TReconCluster.hxx>

namespace CP {
    /// This finds the minimum distance between two clusters, but ignores
    /// outlier hits.  It loops through the hits in the clusters and looks to
    /// see how close the clusters actually get, but ignores up to nDist hits
    /// to calculate the distance.
    double ClusterDistance(const CP::TReconBase& a, 
                           const CP::TReconBase& b,
                           const int nDist = 5);

    /// This finds the minimum distance between hits in the two clusters.
    /// loops through the hits in the clusters and looks to see how close the
    /// clusters actually get.  The cluster also provides the "elliptical"
    /// moments of the cluster that can be used to estimate how close the
    /// clusters get.  That method is faster than ClusterDistance, but is
    /// susceptible to single hit outliers.
    double MinimumClusterDistance(const CP::TReconBase& a, 
                                  const CP::TReconBase& b);

    /// Return a lower bound on the minimum distance between two clusters.
    /// This fast and can be used to check if the clusters are close enough to
    /// bother running the more time consuming ClusterDistance routines.  This
    /// works by placing an XYZ bounding box around each of the clusters and
    /// finds the distance between the boxes.
    double ClusterVicinity(const CP::TReconCluster& a, 
                           const CP::TReconCluster& b);
};

#endif
