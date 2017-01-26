#ifndef TDestroyShortTracks_hxx_seen
#define TDestroyShortTracks_hxx_seen

#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>
#include <TReconTrack.hxx>


namespace CP {
    class TDestroyShortTracks;
    class TReconCluster;
};

/// This takes a algorithm result as a TReconObjectContainer and find
/// and destroy reconstructed tracks with length less than fminLength
/// (now 15mm) into 3D Hits.The 3D hits are returned in allHits
/// container for later use.

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
  
  // Reconstracted tracks with length less than fminLength will not pass this algorihm. Defined in parameters as captRecon.destroyShortTracks.minLength
  
  double fminLength;
  
};

#endif



