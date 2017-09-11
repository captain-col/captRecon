#ifndef TDisassociateHits_hxx_seen
#define TDisassociateHits_hxx_seen

#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>
#include <TReconTrack.hxx>

namespace CP {
    class TDisassociateHits;
    class TReconCluster;
};

/// Save the "big" tracks and then cluster the remaining hits to summarize
/// un-tracked charge distributions.  This takes an algorithm result with a
/// TReconObjectContainer of tracks, showers and clusters.  The goal here is
/// to find groupings of hits that are not consistent with the track
/// hypothesis used in the rest of the reconctruction and then to summarin
/// them as either clusters or showers.  selected objects into their
/// individual hits.  The individual hits are then density clustered for
/// further processing.  Any object that is not broken into its constiuent
/// hits is passed through to the output.
class CP::TDisassociateHits
    : public CP::TAlgorithm {
public:
    TDisassociateHits();
    virtual ~TDisassociateHits();

    /// Apply the algorithm.
    CP::THandle<CP::TAlgorithmResult> 
    Process(const CP::TAlgorithmResult& input,
            const CP::TAlgorithmResult& input1 = CP::TAlgorithmResult::Empty,
            const CP::TAlgorithmResult& input2 = CP::TAlgorithmResult::Empty);

private:
};
#endif
