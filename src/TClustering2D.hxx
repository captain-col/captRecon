#ifndef TClustering2D_hxx_seen
#define TClustering2D_hxx_seen

#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>

namespace CP {
    class TClustering2D;
    class TReconCluster;
};


/// This takes a algorithm result as a TReconObjectContainer and gets all pseudo_3D hits, formed in THitTransfer algorithm. Then it cluster them using DBScan algorithm in to 2D clusters for each plane (X,U,V).

class CP::TClustering2D : public CP::TAlgorithm {
public:
  
    TClustering2D();
  
    virtual ~TClustering2D();

    /// Apply the algorithm.
    CP::THandle<CP::TAlgorithmResult> 
    Process(const CP::TAlgorithmResult& input,
            const CP::TAlgorithmResult& input1 = CP::TAlgorithmResult::Empty,
            const CP::TAlgorithmResult& input2 = CP::TAlgorithmResult::Empty);

  

private:

  unsigned int fminPoints;

  double fmaxDist; 
  
};



#endif
