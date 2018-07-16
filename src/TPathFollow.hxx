#ifndef TPathFollow_hxx_seen
#define TPathFollow_hxx_seen

#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>

namespace CP {
    class TPathFollow;
};


/// Path following algorithm 

class CP::TPathFollow : public CP::TAlgorithm {
public:
  
    TPathFollow();
  
    virtual ~TPathFollow();

    /// Apply the algorithm.
    CP::THandle<CP::TAlgorithmResult> 
    Process(const CP::TAlgorithmResult& input,
            const CP::TAlgorithmResult& input1 = CP::TAlgorithmResult::Empty,
            const CP::TAlgorithmResult& input2 = CP::TAlgorithmResult::Empty);

  

private:

  

};



#endif
