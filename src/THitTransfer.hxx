#ifndef THitTransfer_hxx_seen
#define THitTransfer_hxx_seen
#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>

namespace CP {
    class THitTransfer;
};

/// Modifide Cluster3D pakage, so it takes all 2D hits for each plane and convert tham in to pseudo 3D hits
class CP::THitTransfer
    : public CP::TAlgorithm {
public:
    THitTransfer();
    virtual ~THitTransfer();


    CP::THandle<CP::TAlgorithmResult> 
    Process(const CP::TAlgorithmResult& input,
            const CP::TAlgorithmResult& input1 = CP::TAlgorithmResult::Empty,
            const CP::TAlgorithmResult& input2 = CP::TAlgorithmResult::Empty);


    /// Determine the time zero for the event based on the pmt and wire hits.
    /// The time zero is essentially just time time of the first PMT time
    /// cluster.
    double TimeZero(const CP::THitSelection& pmts,
                    const CP::THitSelection& wires);

private:
  ///Parameter in THitTransfer to reject hits with charge less then this ammoun
  ///They are still get saved, but are not being used in reconstruction at all.
  double fMinCharge;
  
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
 
        
    /// The (up to) three wire hits and make a 3D TWritableReconHit handle.
    /// The hits need to be from different plans (i.e. the wires cannot be
    /// parallel).  This will return an empty handle if there is a problem
    /// constructing the hit.
    bool MakeHit(CP::THitSelection& writableHits,
                 const TVector3& hitPosition,
                 double t0,
                 const CP::THandle<CP::THit>& hit1) const;



};
#endif
