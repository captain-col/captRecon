#ifndef TTrackFitBase_hxx_seen
#define TTrackFitBase_hxx_seen

#include <TReconTrack.hxx>
#include <THandle.hxx>

namespace CP {
    class TTrackFitBase;
};

/// A base class for track fitting. Classes derived from this class implement
/// a specific track fitting algorithm.  There is also a "main" TTrackFit class
/// which serves as a switch yard to determine the best fitter for each type of
/// track.  Notice that this is fitting TReconTrack objects, not TReconPID
/// objects.  TReconPID objects must be fit with a different class.
///
/// This fits a skeleton of a track.  The track is expected to have nodes
/// constructed with a CP::TTrackState and an object derived from
/// CP::TReconBase.  The nodes must be in order from one end of the track to
/// the other.  The input track is expected to be modified by the fitter so
/// that the result handle will be equal to the input handle.  However, this
/// is not guarranteed.  The result track may be a different object than the
/// input track.  If the fit fails, this returns a NULL handle.
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
/// std::auto_ptr<CP::TTrackFitBase> trackFit(new TTrackFit);
/// THandle<TReconTrack> fittedTrack = trackFit->Apply(inputTrack);
/// \endcode
///
/// \warning The input track is expected to be modified by the fitter so that
/// the result handle will be equal to the input handle.  However, this is not
/// guarranteed.  The result track may be a different object than the input
/// track.
class CP::TTrackFitBase {
public:
    virtual ~TTrackFitBase() {}

    /// Fit the skeleton of a track.  This returns NULL if the fit fails.  The
    /// input track is modified during the fit.  See
    /// CP::TTrackFitBase::Apply() for details.
    CP::THandle<CP::TReconTrack>
    operator ()(CP::THandle<CP::TReconTrack>& input) {return Apply(input);}

    /// Fit the skeleton of a track.  The track is expected to have nodes
    /// constructed with a CP::TTrackState and an object derived from
    /// CP::TReconBase.  The nodes must be in order from one end of the track
    /// to the other.  The input track is expected to be modified by the
    /// fitter so that the result handle will be equal to the input handle.
    /// However, this is not guarranteed.  The result track may be a different
    /// object than the input track.  If the fit fails, this returns a NULL
    /// handle.
    virtual CP::THandle<CP::TReconTrack>
    Apply(CP::THandle<CP::TReconTrack>& input) = 0;
};

#endif
