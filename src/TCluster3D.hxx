#ifndef TCluster3D_hxx_seen
#define TCluster3D_hxx_seen
#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>

namespace CP {
    class TCluster3D;
};

/// Take a hit selection (in the form of a CP::TAlgorithmResult), and group
/// the 2D wire hits into 3D hits.  The hits are then formed into a
/// TReconCluster.  If there is a PMT hit selection present, the hits will be
/// placed at a Z associated with the T0 from the PMTs.  If there is no PMT,
/// the clustered hits are all placed at Z=0, with the time set to when the
/// charge passed that plane.
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
    ///
    /// The second and third parameters are ignored by this algorithm.
    CP::THandle<CP::TAlgorithmResult> 
    Process(const CP::TAlgorithmResult& input,
            const CP::TAlgorithmResult& input1 = CP::TAlgorithmResult::Empty,
            const CP::TAlgorithmResult& input2 = CP::TAlgorithmResult::Empty);

    /// Determine the XY crossing point for two wires.  This will throw and
    /// exception if the wires are parallel (e.g. two X wires).  The
    /// z-position of the position is always set to zero.
    TVector3 PositionXY(const CP::THandle<CP::THit>& hit1,
                        const CP::THandle<CP::THit>& hit2);

    /// Determine the time zero for the event based on the pmt and wire hits.
    /// The time zero is essentially just time time of the first PMT time
    /// cluster.
    double TimeZero(const CP::THitSelection& pmts,
                    const CP::THitSelection& wires);

private:

    /// The maximum drift distance in the TPC.  This is set using the
    /// captRecon.Cluster3D.maxDrift
    double fMaxDrift;

    /// The time around the central time of the X hit where a V or U hit is
    /// considered to overlap.  This is in units of the time RMS.
    double fXSeparation;

    /// The time around the central time of the V hit where a X or U hit is
    /// considered to overlap.  This is in units of the time RMS.
    double fVSeparation;

    /// The time around the central time of the U hit where a V or X hit is
    /// considered to overlap.  This is in units of the time RMS.
    double fUSeparation;

    /// Hits that are closer than this in time are considered to overlap.
    /// This protects against unreasonably small RMSs for the hit times.
    double fMinSeparation;
};
#endif
