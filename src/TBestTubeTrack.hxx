#ifndef TBestTubeTrack_hxx_seen
#define TBestTubeTrack_hxx_seen
#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>

namespace CP {
    class TBestTubeTrack;
};

/// This takes a algorithm result with clusters sliced by Z (e.g. from
/// TClusterSlice), and turns them into tracks.  The tracks are found by
/// running a Best Tube algorithm looking at random pairs of clusters in 3D.
/// The track is then assumed to go through the center of each cluster.  This
/// does pattern recognition, and then does a first approximation fit to a
/// track.
class CP::TBestTubeTrack
    : public CP::TAlgorithm {
public:
    TBestTubeTrack();
    virtual ~TBestTubeTrack();

    /// Apply the algorithm.
    CP::THandle<CP::TAlgorithmResult> 
    Process(const CP::TAlgorithmResult& input,
            const CP::TAlgorithmResult& input1 = CP::TAlgorithmResult::Empty,
            const CP::TAlgorithmResult& input2 = CP::TAlgorithmResult::Empty);

private:
    /// Take an object container of clusters created by slicing the hits into
    /// "z-slices" (the output of MakeSlices) and split them up into tracks.
    /// The output is a object container that may contain one or more tracks.
    /// This can take a handle to the result being created that can be filled
    /// with debugging information about the track finding will be filled into
    /// the output result.
    CP::THandle<CP::TReconObjectContainer> 
    MakeTracks(const CP::TReconObjectContainer& input, 
               CP::THandle<CP::TAlgorithmResult> thisResult);

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

};
#endif
