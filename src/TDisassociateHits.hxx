#ifndef TDisassociateHits_hxx_seen
#define TDisassociateHits_hxx_seen

#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>
#include <TReconTrack.hxx>

namespace CP {
    class TDisassociateHits;
    class TReconCluster;
};

/// This takes a algorithm result with a TReconObjectContainer of with objects
/// and disassociates selected objects into their individual hits.  The
/// individual hits are then density clustered for further processing.  Any
/// object that is not broken into its constiuent hits is passed through to
/// the output.
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
