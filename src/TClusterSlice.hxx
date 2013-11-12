#ifndef TClusterSlice_hxx_seen
#define TClusterSlice_hxx_seen
#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>

namespace CP {
    class TClusterSlice;
};

/// This takes a algorithm result with 3D hits and then splits the hits into
/// clusters in "Z".  This means that 3D hits which are part of the same XY
/// plane confusion region are all in the same cluster.  This simplifies the
/// infirmation that needs to be looked at for tracking.
class CP::TClusterSlice
    : public CP::TAlgorithm {
public:
    TClusterSlice();
    virtual ~TClusterSlice();

    /// Apply the algorithm.
    CP::THandle<CP::TAlgorithmResult> 
    Process(const CP::TAlgorithmResult& input,
            const CP::TAlgorithmResult& input1 = CP::TAlgorithmResult::Empty,
            const CP::TAlgorithmResult& input2 = CP::TAlgorithmResult::Empty);

private:
    /// Take hits and break them into clusters based on the "Z Confusion
    /// Distance".  The slices along the Z axis are controlled by the
    /// captRecon.clusterSlice.clusterStep parameter.  The output is an object
    /// container of 3D aclusters.
    CP::THandle<CP::TReconObjectContainer> 
    MakeSlices(CP::THandle<CP::THitSelection> input);

    /// The minimum number of neighbors within the maximum distance of the
    /// current point to consider the current point to be in a high density
    /// region.
    int fMinPoints;

    /// The radius over which points are counted to determine if a point is in
    /// a high density region.
    int fMaxDist;

    /// The minimum number of hits before a cluster is broken up.  This keeps
    /// "micro" clusters from being split up since they will be best handled
    /// using the cluster directly.
    unsigned int fMinHits;

    /// The splitting distance in Z.
    double fClusterStep;
};
#endif
