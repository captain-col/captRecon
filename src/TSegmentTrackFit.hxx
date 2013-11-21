#ifndef TSegmentTrackFit_hxx_seen
#define TSegmentTrackFit_hxx_seen

#include "TTrackFitBase.hxx"

#include <THitSelection.hxx>
#include <TReconTrack.hxx>
#include <THandle.hxx>

namespace CP {
    class TSegmentTrackFit;
};

/// Fit a track by forming segments between each of the track node objects.
/// The resulting track is a sequence of straight line segments passing
/// through the center of the track node objects.  This fit is mostly useful
/// where the TReconTrack object is being used to hold a sequence of clusters
/// where the order matters, but it is being used as an intermediate result in
/// pattern recognition.  The direction for each node is the average of the
/// vectors to the objects on either side. This takes an input track where all
/// of the nodes are filled with objects, and the nodes are in the correct
/// order.  See the CP::TTrackFitBase class for more detailed API
/// documentation.
class CP::TSegmentTrackFit : public CP::TTrackFitBase {
public:
    TSegmentTrackFit();
    virtual ~TSegmentTrackFit();

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
