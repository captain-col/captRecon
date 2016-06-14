#ifndef TClusterMerge_hxx_seen
#define TClusterMerge_hxx_seen
#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>
#include <TReconCluster.hxx>
#include "TPositionNeighbors.hxx"

namespace CP {
    class TClusterMerge;
};

/// Take a collection of clusters as produced by TClusterSlice or similar
/// algorithm (in the form of a CP::TAlgorithmResult), and merge groups that
/// have been split along the flat axis.  This algorithm is intended to follow
/// an algorithm that has clustered the hits into time slices, and to undo
/// splitting that happens for horizontal tracks, and to resplit broad
/// clusters of hits along the preferred axis.
class CP::TClusterMerge
    : public CP::TAlgorithm {
    typedef CP::TPositionNeighbors< CP::THandle<CP::THit> > Neighbors;
public:
    TClusterMerge();
    virtual ~TClusterMerge();

    /// Take a TAlgorithmResult handle containing a container of clusters and
    /// optimize the clustering.  It produces a TReconObject container of
    /// clusters that are ready for track formation.  This usually follows
    /// TClusterSlice.
    ///
    /// Produced object containers
    ///
    /// * final : The main output container.
    ///
    /// * merged : A container of merged clusters before they are split up to
    ///            help form tracks.
    ///
    /// The second and third parameters are ignored by this algorithm.
    CP::THandle<CP::TAlgorithmResult> 
    Process(const CP::TAlgorithmResult& input,
            const CP::TAlgorithmResult& input1 = CP::TAlgorithmResult::Empty,
            const CP::TAlgorithmResult& input2 = CP::TAlgorithmResult::Empty);

private:
    /// Combine two clusters into a new cluster.  The original clusters are not
    /// modified.
    CP::THandle<CP::TReconCluster> CombineClusters(
        const CP::THandle<CP::TReconCluster>& cluster1,
        const CP::THandle<CP::TReconCluster>& cluster2);

    /// Check if two clusters are overlapping along the direction of a track.
    /// This is looking for clusters that have been sliced at a high angle to
    /// the track direction (e.g. cluster slices for a horizontal track)
    bool OverlappingClusters(
        CP::TClusterMerge::Neighbors& neighbors,
        const CP::THandle<CP::TReconCluster>& cluster1,
        const CP::THandle<CP::TReconCluster>& cluster2);

    /// The minimum size for a cluster to be checked if it should merge with a
    /// neighbor.
    int fMinimumClusterSize;

    /// The width of a 3D hit along the time axis (i.e. Z) in units of RMS.
    /// This is set in terms of the RMS times a multiplier, with a value of
    /// 1.0 meaning the hit is assumed to have an extent of 1.0 times the RMS.
    double fClusterSizeFactor;

    /// If either of two clusters overlap with a neighbor by more than this
    /// amount, then the two clusters are taken as overlapping.
    double fClusterOverlap;

    /// Set the maximum width as a fraction of the cluster width for a cluster
    /// to be split after merging.  After clusters have been merged, they can
    /// be re-split along the cluster length.
    double fMaximumWidthFraction;

    /// Set the minimum length for a cluster to be resplit along the long
    /// axis.  Clusters with a long axis shorter than this length will remain
    /// as clusters (i.e. a small stand alone blip of energy), or be
    /// incorporated into another track/shower.
    double fMinimumLength;

    /// The extent of the cluster to be used when splitting a merged cluster
    /// into "optimal" sub-clusters.  This is the extent in the XY plane
    /// (i.e. the "plane of confusion").
    double fClusterExtent;

    /// The distance metric to be used along the Z direction (i.e. the time
    /// axis).  The metric along the X and Y axis will be 1.0.
    double fTimeMetric;
    
};
#endif
