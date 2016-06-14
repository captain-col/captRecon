#ifndef TClusterSlice_hxx_seen
#define TClusterSlice_hxx_seen
#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>

namespace CP {
    class TClusterSlice;
};

/// This takes a algorithm result with 3D hits and then splits the hits into
/// clusters in "Z".  This means that 3D hits which are part of the same XY
/// plane confusion region are all in the same cluster.  The hits in the "Z"
/// slice are then clustered using the dbscan algorithm.  This simplifies the
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

    /// The minimum number of hits in the input in the input object for the
    /// hits to be sliced into clusters.  This keeps "micro" events from
    /// being split up since they will be best handled using the cluster
    /// directly.  This is set using captRecon.clusterSlice.minHits.
    unsigned int fMinHits;

    /// The minimum number hits in any given slice and the extent of the slice
    /// will continue to grow until this criteria is met.  This is set using
    /// captRecon.clusterSlice.minPoints
    int fMinPoints;

    /// A minimum distance to be made into a slice.  This prevents
    /// ridiculously small slices in a very big event.  It is set using
    /// captRecon.clusterSlice.minStep.
    double fMinStep;

    /// The target thickness in Z for the slices.  The actual thickness is
    /// then calculated to have an integer number of slices.  This is set
    /// using captRecon.clusterSlice.clusterStep.
    double fClusterStep;

    /// The splitting distance in X/Y.  This is used to define the DBScan
    /// metric used to separate a slice into localized clusters.  This is set
    /// captRecon.clusterSlice.clusterExtent.
    double fClusterExtent;

    /// The scaling factor for the cluster growth.
    double fClusterGrowth;
    
    /// The required amount of charge to be in any cluster.  This prevents
    /// really small clusters from being formed.  This is not parameterized,
    /// since it sets a very low floor on the cluster size.
    double fClusterCharge;
};
#endif
