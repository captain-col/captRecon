#ifndef TTracking2D_hxx_seen
#define TTracking2D_hxx_seen

#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>

namespace CP {
    class TTracking2D;
};


/// This takes a algorithm result as a TReconObjectContainer and gets all formed clusters and convert tham in to tracks. It should be a cut of length/width of the cluster ratio, but apparantly it does not work, so it is skipped.

class CP::TTracking2D : public CP::TAlgorithm {
public:
  
    TTracking2D();
  
    virtual ~TTracking2D();

    /// Apply the algorithm.
    CP::THandle<CP::TAlgorithmResult> 
    Process(const CP::TAlgorithmResult& input,
            const CP::TAlgorithmResult& input1 = CP::TAlgorithmResult::Empty,
            const CP::TAlgorithmResult& input2 = CP::TAlgorithmResult::Empty);

  

private:

  double fXZRatio;

};



#endif
