#ifndef TMinimalSpanningTrack_hxx_seen
#define TMinimalSpanningTrack_hxx_seen
#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>

namespace CP {
    class TMinimalSpanningTrack;
};

/// This takes a algorithm result with clusters (e.g. from TClusterSlice), and
/// turns them into tracks.  The tracks are found by finding a minimal
/// spanning tree for all of the clusters based on the distance between the
/// clusters.  The spanning tree is then turned into a (possibly) large set of
/// tracks that are intended for later processing, plus any clusters that
/// don't get added to a track.  The tracks are fit using the
/// TSegmentTrackFit.
///
/// \warning The result TReconObjectContainer (in the TAlgorithmResult) will
/// probably contain both tracks and clusters.
class CP::TMinimalSpanningTrack
    : public CP::TAlgorithm {
public:
    TMinimalSpanningTrack();
    virtual ~TMinimalSpanningTrack();

    /// Apply the algorithm.
    CP::THandle<CP::TAlgorithmResult> 
    Process(const CP::TAlgorithmResult& input,
            const CP::TAlgorithmResult& input1 = CP::TAlgorithmResult::Empty,
            const CP::TAlgorithmResult& input2 = CP::TAlgorithmResult::Empty);

private:
};
#endif
