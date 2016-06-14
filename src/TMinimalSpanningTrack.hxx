#ifndef TMinimalSpanningTrack_hxx_seen
#define TMinimalSpanningTrack_hxx_seen
#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>

namespace CP {
    class TMinimalSpanningTrack;
};

/// This takes a algorithm result with clusters (e.g. from TClusterSlice or
/// TClusterMerge), and turns them into track candidates.  The candidates are
/// found by forming a minimal spanning tree for all of the clusters based on
/// the distance between the clusters.  The spanning tree is then turned into
/// a (possibly) large set of tracks that are intended for later processing,
/// plus any clusters that don't get added to a track.  The tracks are fit
/// using the TSegmentTrackFit.
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

    /// The maximum distance between cluster centroids for an edge to be
    /// included in the graph.  This is an optimization to help the algorithm
    /// beat the n-squared problem before the hit distance cut is applied.  
    double fDistCut;

    /// The maximum distance between clusters based on the closest hits.  This
    /// is a "slightly" fuzzy cut since the distance is based on
    /// ClusterVicinity and not the actual closest approach.  Any clusters
    /// that are closer than this will have an edge placed into the graph.
    /// This is an optimization to help short circuit the algorithm and reduce
    /// the total number of edges.
    double fHitDistCut;
};
#endif
