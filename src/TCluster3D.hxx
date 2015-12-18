#ifndef TCluster3D_hxx_seen
#define TCluster3D_hxx_seen
#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>

namespace CP {
    class TCluster3D;
};

/// Take a hit selection (in the form of a CP::TAlgorithmResult), and group
/// the 2D wire hits into 3D hits that are returned as a THitSelection (in the
/// output TAlgorithm Result.  The 3D hits form a hit cloud representing the
/// best knowledge of the spacial distribution of drifting electrons.  If
/// there is an input PMT hit selection present, the hits will be placed at a
/// Z associated with the T0 from the PMTs.  If there is no PMT, T0 will be
/// taken to be zero.
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

    /// Determine the XY crossing point for two wires.  This will throw and
    /// exception if the wires are parallel (e.g. two X wires).  The
    /// z-position of the position is always set to zero.
    TVector3 PositionXY(const CP::THandle<CP::THit>& hit1,
                        const CP::THandle<CP::THit>& hit2);

    /// Determine if three hits overlap at a single point.  If the wires
    /// overlap, the position is returned by reference.  The z-position of the
    /// result is always set to zero.
    bool OverlapXY(const CP::THandle<CP::THit>& hit1,
                   const CP::THandle<CP::THit>& hit2,
                   const CP::THandle<CP::THit>& hit3,
                   TVector3& position);

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

    /// A 3D hit is formed from 3 wire hits which are observing the same
    /// charge distribution.  Since each of the three wire hits may have
    /// contributions from multiple "electron clouds", the most charge that a
    /// single 3D hit can have is the minimum charge on any single wire (it
    /// might have less when the charge is properly shared between 3D hits).
    /// This means that the charge on the smallest wire hit is a fairly good
    /// estimator for the expected charge in the 3D hit.  A better estimate of
    /// the 3D hit charge is look at the charge in the samples which overlap
    /// in time.  Hits where the overlapping charge is less than
    /// fMinimumOverlap times the expected charge are rejected.  This rejects
    /// wires that are measureing drift clouds at different times.
    double fMinimumOverlap;

    /// Hits that are closer than this in time are considered to overlap.
    /// This protects against unreasonably small RMSs for the hit times.
    double fMinSeparation;

    /// Set the maximum time spread for a reconstructed 3D hits.  This keeps a
    /// 3D hit from being long and skinny along the time axis (i.e. Z).
    double fMaximumSpread;

    /// The amount of energy per measured ionization electron.  This is mostly
    /// used for pretty output and diagnostics.  The hits have their charges
    /// in electrons
    double fEnergyPerCharge;

    /// The amount of time per digitizer sample.  This is the same for all
    /// hits (and is nominally 500*ns).
    double fDigitStep;
    
    /// Take the "size" of each hit, and the minimum separation (basically the
    /// time bin size), and return an overlap time.  If hits are within the
    /// overlap time of each other, then they are considered to overlap.
    double OverlapTime(double rms1, double rms2, double step) const;

    /// Return true if two hits overlap.
    bool OverlappingHits(CP::THandle<CP::THit> h1,
                         CP::THandle<CP::THit> h2) const;
        
    /// The (up to) three wire hits and make a 3D TWritableReconHit handle.
    /// The hits need to be from different plans (i.e. the wires cannot be
    /// parallel).  This will return an empty handle if there is a problem
    /// constructing the hit.
    bool MakeHit(CP::THitSelection& writableHits,
                 const TVector3& hitPosition,
                 double t0,
                 const CP::THandle<CP::THit>& hit1,
                 const CP::THandle<CP::THit>& hit2,
                 const CP::THandle<CP::THit>& hit3) const;

    /// Find the fractional overlap of the hit with the constituent.
    double FindOverlap(const CP::THandle<CP::THit>& hit,
                       const CP::THandle<CP::THit>& constituent) const;

};
#endif
