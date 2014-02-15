#ifndef TBootstrapTrackFit_hxx_seen
#define TBootstrapTrackFit_hxx_seen

#include "TTrackFitBase.hxx"
#include "ECaptRecon.hxx"

#include <THitSelection.hxx>
#include <TReconTrack.hxx>
#include <THandle.hxx>
#include <ECore.hxx>

namespace CP {

    /// A base fitting object.
    EXCEPTION(EBootstrapTrackFit,ECaptRecon);

    /// A base fitting object.
    EXCEPTION(EBootstrapMissingNodeObject,EBootstrapTrackFit);

    class TBootstrapTrackFit;
};

/// Fit a track using the BFL Bootstrap particle filter algorithm.  This takes
/// an input track where all of the nodes are filled with objects, and the
/// nodes are in the correct order.  See the CP::TTrackFitBase class for more
/// detailed API documentation.  The constructor takes a single argument which
/// is the number of trials for each iteration.  It has a default value.
class CP::TBootstrapTrackFit : public CP::TTrackFitBase {
public:

    /// The argument is the number of trials at each iteration.
    explicit TBootstrapTrackFit(int trials=500);
    virtual ~TBootstrapTrackFit();

    /// Fit the skeleton of a track using a bootstrap particle filter.  This
    /// is strictly speaking a "smoother" since it will using forward/backward
    /// smoothing.
    virtual CP::THandle<CP::TReconTrack>
    Apply(CP::THandle<CP::TReconTrack>& input);

private:
    /// Estimate the range of the track.
    double EstimateRange(const CP::TReconNodeContainer& nodes);

    /// Estimate the momentum of the track (needed to estimate multiple
    /// scattering.   This assumes a muon track so it's a really crude estimate.
    double EstimateMomentum(const CP::TReconNodeContainer& nodes);

    /// The number of trials for each iteration of the fitter.
    int fTrials;
};

#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// End:
