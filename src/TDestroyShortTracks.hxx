#ifndef TDestroyShortTracks_hxx_seen
#define TDestroyShortTracks_hxx_seen

#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>
#include <TReconTrack.hxx>


namespace CP {
    class TDestroyShortTracks;
    class TReconCluster;
};

/// This takes a algorithm result with a TReconObjectContainer of with objects
/// and disassociates selected objects into their individual hits.  The
/// individual hits are then density clustered for further processing.  Any
/// object that is not broken into its constiuent hits is passed through to
/// the output.
class CP::TDestroyShortTracks
    : public CP::TAlgorithm {
public:
    TDestroyShortTracks();
    virtual ~TDestroyShortTracks();

    /// Apply the algorithm.
    CP::THandle<CP::TAlgorithmResult> 
    Process(const CP::TAlgorithmResult& input,
            const CP::TAlgorithmResult& input1 = CP::TAlgorithmResult::Empty,
            const CP::TAlgorithmResult& input2 = CP::TAlgorithmResult::Empty);

private:
};
#endif



