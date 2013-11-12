#ifndef TClusterTrackFit_hxx_seen
#define TClusterTrackFit_hxx_seen

#include "TTrackFitBase.hxx"

#include <THitSelection.hxx>
#include <TReconTrack.hxx>
#include <THandle.hxx>

namespace CP {
    class TClusterTrackFit;
};

/// Fit a track by forming a cluster from the track node objects.  This takes
/// an input track where all of the nodes are filled with objects, and the
/// nodes are in the correct order.  See the CP::TTrackFitBase class for more
/// detailed API documentation.
class CP::TClusterTrackFit : public CP::TTrackFitBase {
public:
    TClusterTrackFit();
    virtual ~TClusterTrackFit();

    /// Fit the skeleton of a track using the moments of a cluster formed from
    /// the track node objects.
    virtual CP::THandle<CP::TReconTrack>
    Apply(CP::THandle<CP::TReconTrack>& input);
};

#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// End:
