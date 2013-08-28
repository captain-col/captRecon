#ifndef TTrackFit_hxx_seen
#define TTrackFit_hxx_seen

#include "TTrackFitBase.hxx"

#include <THitSelection.hxx>
#include <TReconTrack.hxx>
#include <THandle.hxx>

namespace CP {
    class TTrackFit;
};

/// This is a wrapper around other track fitting classes (all derived from
/// TTrackFitBase that chooses the correct fitter to be applied.  The
/// TTrackFit class is the "main" class serving as a switch yard to determine
/// the best fitter for each type of track.  Notice that this is fitting
/// TReconTrack objects, not TReconPID objects.  TReconPID objects must be fit
/// with a different class.
///
/// Most code should be using TTrackFit which will choose the best fitter to
/// use in each circumstance.  How to use these fitting classes:
/// \code
/// TTrackFit trackFit;
/// THandle<TReconTrack> fittedTrack = trackFit(inputTrack);
/// if (fittedTrack) std::cout << "fit was successful" << std::endl;
/// if (!fittedTrack) std::cout << "fit failed" << std::endl;
/// \endcode
/// If the fit fails then the returned handle will be empty.  The track fit
/// can also be applied using the Apply method which will be more convenient
/// when the fitter is referenced by a pointer.
/// \code
/// std::auto_ptr<CP::TTrackFit> trackFit(new TTrackFit);
/// THandle<TReconTrack> fittedTrack = trackFit->Apply(inputTrack);
/// \endcode
///
/// \warning The input track is expected to be modified by the fitter so that
/// the result handle will be equal to the input handle.  However, this is not
/// guarranteed.  The result track may be a different object than the input
/// track.
class CP::TTrackFit : public CP::TTrackFitBase {
public:
    TTrackFit();
    virtual ~TTrackFit();

    /// Fit the skeleton of a track.  The track is expected to have nodes
    /// constructed with a CP::TTrackState and an object derived from
    /// CP::TReconBase.  The nodes must be in order from one end of the track
    /// to the other.  The input track is expected to be modified by the
    /// fitter so that the result handle will be equal to the input handle.
    /// However, this is not guarranteed.  The result track may be a different
    /// object than the input track.  If the fit fails, this returns a NULL
    /// handle.
    virtual CP::THandle<CP::TReconTrack>
    Apply(CP::THandle<CP::TReconTrack>& input);
};

#endif
