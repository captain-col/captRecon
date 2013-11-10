#ifndef TClusterSlice_hxx_seen
#define TClusterSlice_hxx_seen
#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>

namespace CP {
    class TClusterSlice;
};

/// This takes a algorithm result with clusters of "simply connected" hits
/// (meaning hits that are connected by direct neighbors), and then splits
/// each cluster into clusters in "Z".  The resulting clusters are then turned
/// into tracks.  This means that 3D hits which are part of the same XY plane
/// confusion region are all in the same cluster.  The track is then assumed
/// to go through the center of each cluster.  This is intended as a first
/// approximation to a track.
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
    /// Take a single cluster of simply connected hits and turn it into a
    /// track.  If the hits are not long enough to turn into a track, then
    /// the original cluster is returned.  If the input object is not a
    /// cluster, then the original object is returned.
    CP::THandle<CP::TReconBase> 
    ClusterToTrack(CP::THandle<CP::TReconBase> input);

    /// Take hits and break them into clusters based on the "Z Confusion
    /// Distance".  The slices along the Z axis are controlled by the
    /// captRecon.clusterSlice.clusterStep parameter.  The output is an object
    /// container of 3D aclusters.
    CP::THandle<CP::TReconObjectContainer> 
    MakeSlices(CP::THandle<CP::THitSelection> input);

    /// Take an object container of clusters created by slicing a cluster of
    /// simply connected hits into "z-slices" (the output of BreakCluster) and
    /// split them up into tracks.  The output is a object container that may
    /// contain one or more tracks.  It can have more than one track when the
    /// simply connected cluster is actually two tracks connected to a vertex
    /// (or some other similar topology).  The seeds container will be filled
    /// with the seeds that were found.
    CP::THandle<CP::TReconObjectContainer> 
    MakeTracks(const CP::TReconObjectContainer& input,
               CP::TReconObjectContainer& seeds);

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
