#ifndef TClusterUnusedHits_hxx_seen
#define TClusterUnusedHits_hxx_seen

#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>
#include <TReconTrack.hxx>

namespace CP {
    class TClusterUnusedHits;
    class TReconCluster;
};


/// This takes a algorithm result as a TReconObjectContainer and all hits
/// in the event(THitSelection). The algorifm finds 3D hits, which were not
/// used in andy object from TReconObjectContainer and cluter tham by position.
/// Formed clusters are saved in oucome as well as all objects from
/// input TReconObjectContainer.

class CP::TClusterUnusedHits
    : public CP::TAlgorithm {
public:
    TClusterUnusedHits();
    virtual ~TClusterUnusedHits();

    /// Apply the algorithm.
    CP::THandle<CP::TAlgorithmResult> 
    Process(const CP::TAlgorithmResult& input,
            const CP::TAlgorithmResult& input1,
            const CP::TAlgorithmResult& input2 = CP::TAlgorithmResult::Empty);

private:

  unsigned int fminPoints;

  double fmaxDist; 
  
};



#endif
