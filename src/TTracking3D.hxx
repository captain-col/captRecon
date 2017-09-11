#ifndef TTracking3D_hxx_seen
#define TTracking3D_hxx_seen

#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>

namespace CP {
    class TTracking3D;
};


/// This takes a algorithm result as a TReconObjectContainer gets all 2D tracks, finds 3(or 2 but always with Xplane track) tracks happened in on Zcoodrinate frame. Then the algorithm devide them in hits along the way and combine hits with appropriate Z in 3D clusters wich will serve to create 3D tracks.

class CP::TTracking3D : public CP::TAlgorithm {
public:
  
    TTracking3D();
  
    virtual ~TTracking3D();

    /// Apply the algorithm.
    CP::THandle<CP::TAlgorithmResult> 
    Process(const CP::TAlgorithmResult& input,
            const CP::TAlgorithmResult& input1 = CP::TAlgorithmResult::Empty,
            const CP::TAlgorithmResult& input2 = CP::TAlgorithmResult::Empty);

  

private:


};



#endif
