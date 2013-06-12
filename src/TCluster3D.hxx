#ifndef TCluster3D_hxx_seen
#define TCluster3D_hxx_seen
#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>

namespace CP {
    class TCluster3D;
};

/// Take a hit selection (in the form of a CP::TAlgorithmResult), and group
/// the 2D wire hits into 3D hits.
class CP::TCluster3D
    : public CP::TAlgorithm {
public:
    TCluster3D();
    virtual ~TCluster3D();

    /// Take a TAlgorithmResult handle containing 2D wire hits, and group them
    /// into 3D XYZT hits.  Since a THitSelection can be converted into a
    /// TAlgorithmResult this can also be called with a THandle to a
    /// THitSelection.  The output TAlgorithmResult will contain:
    ///
    ///   * used  -- All of the wire hits that were combined into 3D hits.
    /// 
    ///   * unused -- Any wire hits that were not used in 3D hits.
    ///
    ///   * clustered -- The hit selection containing the 3D TReconHit from
    ///                  this algorithm.  This is the last THitSelection
    ///                  added.
    CP::THandle<CP::TAlgorithmResult> 
    Process(const CP::TAlgorithmResult& input);

    /// Determine the XY crossing point for two wires.  This will throw and
    /// exception if the wires are parralel (e.g. two X wires).
    TVector3 PositionXY(const CP::THandle<CP::THit>& hit1,
                        const CP::THandle<CP::THit>& hit2);
};
#endif
