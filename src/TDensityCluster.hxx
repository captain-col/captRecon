#ifndef TDensityCluster_hxx_seen
#define TDensityCluster_hxx_seen
#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>

namespace CP {
    class TDensityCluster;
};

/// Take a hit selection of 3D hits, and cluster the hits into contiguous
/// clusters using a density cluster algorithm.  The goal is to find simply
/// connected sets of hits that can then be analyzed as tracks (or possibly
/// showers).  The core of this algorithm uses the TTmplDensityCluster
/// implementation of DBSCAN.
class CP::TDensityCluster
    : public CP::TAlgorithm {
public:
    TDensityCluster();
    virtual ~TDensityCluster();

    /// Apply the algorithm.
    CP::THandle<CP::TAlgorithmResult> 
    Process(const CP::TAlgorithmResult& input,
            const CP::TAlgorithmResult& input1 = CP::TAlgorithmResult::Empty,
            const CP::TAlgorithmResult& input2 = CP::TAlgorithmResult::Empty);

private:
    /// The minimum number of neighbors within the maximum distance of the
    /// current point to consider the current point to be in a high density
    /// region.
    int fMinPoints;

    /// The radius over which points are counted to determine if a point is in
    /// a high density region.
    int fMaxDist;
};
#endif
