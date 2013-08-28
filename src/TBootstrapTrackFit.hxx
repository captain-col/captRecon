#ifndef TBootstrapTrackFit_hxx_seen
#define TBootstrapTrackFit_hxx_seen

#include "TTrackFitBase.hxx"

#include <THitSelection.hxx>
#include <TReconTrack.hxx>
#include <THandle.hxx>

namespace CP {
    class TBootstrapTrackFit;
};

/// Fit a track using the BFL Bootstrap particle filter algorithm.  This takes
/// an input track where all of the nodes are filled with objects, and the
/// nodes are in the correct order.  See the CP::TTrackFitBase class for more
/// detailed API documentation.
class CP::TBootstrapTrackFit : public CP::TTrackFitBase {
public:
    TBootstrapTrackFit();
    virtual ~TBootstrapTrackFit();

    /// Fit the skeleton of a track using a bootstrap particle filter.  This
    /// is strictly speaking a "smoother" since it will using forward/backward
    /// smoothing.
    virtual CP::THandle<CP::TReconTrack>
    Apply(CP::THandle<CP::TReconTrack>& input);
};

#endif
