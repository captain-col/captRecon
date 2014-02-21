#ifndef TTrackMassFit_hxx_seen
#define TTrackMassFit_hxx_seen

#include <TReconTrack.hxx>
#include <THandle.hxx>

namespace CP {
    class TTrackMassFit;
};

/// A class to the mass and total energy deposition for a track that has
/// already been fit using TTrackFit.  The input track is expected to be
/// modified by the fitter so that the result handle will be equal to the
/// input handle.  However, this is not guarranteed.  The result track may be
/// a different object than the input track.  If the fit fails, this returns a
/// NULL handle.
///
/// This is how to use these fitting classes:
/// \code
/// TTrackMassFit massFit;
/// THandle<TReconTrack> fittedTrack = massFit(inputTrack);
/// if (fittedTrack) std::cout << "fit was successful" << std::endl;
/// if (!fittedTrack) std::cout << "fit failed" << std::endl;
/// \endcode
/// If the fit fails then the returned handle will be empty.  The track fit
/// can also be applied using the Apply method which will be more convenient
/// when the fitter is referenced by a pointer.
/// \code
/// std::auto_ptr<CP::TTrackMassFit> massFit(new TTrackMassFit);
/// THandle<TReconTrack> fittedTrack = massFit->Apply(inputTrack);
/// \endcode
///
/// \warning The input track is expected to be modified by the fitter so that
/// the result handle will be equal to the input handle.  However, this is not
/// guarranteed.  The result track may be a different object than the input
/// track.
class CP::TTrackMassFit {
public:
    TTrackMassFit();
    virtual ~TTrackMassFit();

    /// Fit the mass and energy deposition for track.  This returns NULL if
    /// the fit fails.  The input track is modified during the fit.  See
    /// CP::TTrackMassFit::Apply() for details.
    CP::THandle<CP::TReconTrack>
    operator ()(CP::THandle<CP::TReconTrack>& input) {return Apply(input);}

    /// Fit the mass and energy deposition for track.  The track is expected
    /// to have already been fit using a class derived from TTrackFitBase
    /// (usually TTrackFit).  The input track is expected to be modified by
    /// the fitter so that the result handle will be equal to the input
    /// handle.  However, this is not guarranteed.  The result track may be a
    /// different object than the input track.  If the fit fails, this returns
    /// a NULL handle. 
    ///
    /// If the track sense can be determined (i.e. which direction it's
    /// going), then this fitter will reverse the track so that start of the
    /// track is at the fitted beginning.  This fitter takes into account
    /// whether the track is fully contained, or exiting.
    virtual CP::THandle<CP::TReconTrack>
    Apply(CP::THandle<CP::TReconTrack>& input);
};

#endif
