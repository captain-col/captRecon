#ifndef TFitting3D_hxx_seen
#define TFitting3D_hxx_seen

#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>

#include <set>

namespace CP {
    class TFitting3D;
};

/// Apply Fitter to 3D tracks
class CP::TFitting3D
    : public CP::TAlgorithm {
public:
    typedef std::set< CP::THandle<CP::THit> > HitSet;

    TFitting3D();
    virtual ~TFitting3D();

    /// Apply the algorithm.
    CP::THandle<CP::TAlgorithmResult> 
    Process(const CP::TAlgorithmResult& input,
            const CP::TAlgorithmResult& input1 = CP::TAlgorithmResult::Empty,
            const CP::TAlgorithmResult& input2 = CP::TAlgorithmResult::Empty);

    
private:

};
#endif
