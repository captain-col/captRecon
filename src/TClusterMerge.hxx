#ifndef TClusterMerge_hxx_seen
#define TClusterMerge_hxx_seen
#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>
#include <TReconCluster.hxx>
#include "TPositionNeighbors.hxx"

namespace CP {
    class TClusterMerge;
};

/// Take a collection of clusters (in the form of a CP::TAlgorithmResult), and
/// merge groups that have been split.  This algorithm is intended to follow
/// an algorithm that has clustered the hits into time slices, and to undo
/// splitting that happens for horizontal tracks, and to resplit broad
/// clusters of hits along the prefered axis.
class CP::TClusterMerge
    : public CP::TAlgorithm {
    typedef CP::TPositionNeighbors< CP::THandle<CP::THit> > Neighbors;
public:
    TClusterMerge();
    virtual ~TClusterMerge();

    /// Take a TAlgorithmResult handle containing a container of clusters and
    /// optimize the clustering.  This usually follows TClusterSlice.
    ///
    ///   * used -- A THitSelection of all of the wire hits that were combined
    ///             into 3D hits.
    /// 
    ///   * unused -- A THitSelection of any wire hits that were not used in
    ///               3D hits.
    ///
    ///   * clustered -- The hit selection containing the 3D TReconHit from
    ///                  this algorithm.  This is the last THitSelection
    ///                  added.
    ///
    /// The second and third parameters are ignored by this algorithm.
    CP::THandle<CP::TAlgorithmResult> 
    Process(const CP::TAlgorithmResult& input,
            const CP::TAlgorithmResult& input1 = CP::TAlgorithmResult::Empty,
            const CP::TAlgorithmResult& input2 = CP::TAlgorithmResult::Empty);

private:
    // Combine two clusters into a new cluster.  The original clusters are not
    // modified.
    CP::THandle<CP::TReconCluster> CombineClusters(
        const CP::THandle<CP::TReconCluster>& cluster1,
        const CP::THandle<CP::TReconCluster>& cluster2);

    // Check if two clusters are overlapping along the direction of a track.
    // This is looking for clusters that have been sliced at a high angle to
    // the track direction (e.g. cluster slices for a horizontal track)
    bool OverlappingClusters(
        CP::TClusterMerge::Neighbors& neighbors,
        const CP::THandle<CP::TReconCluster>& cluster1,
        const CP::THandle<CP::TReconCluster>& cluster2);

};
#endif
