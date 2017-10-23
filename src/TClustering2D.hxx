#ifndef TClustering2D_hxx_seen
#define TClustering2D_hxx_seen

#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>

namespace CP {
    class TClustering2D;
    class TReconCluster;
};


/// This takes a algorithm result as a TReconObjectContainer and gets all pseudo_3D hits, formed in THitTransfer algorithm.
/// Then For each plane separately the algorithm applies position based clustering(parameters: minpoints=3,masdist=25mm)<- it is done to take away "noise". Then for all hits that where clusterd the HoughTransform where applyed.
///Hough produce Line that might be splitted in to several lines based on gaps by applying DBScan(parameters: minpoints=2,masdist=60mm)

class CP::TClustering2D : public CP::TAlgorithm {
public:
  
    TClustering2D();
  
    virtual ~TClustering2D();

    /// Apply the algorithm.
    CP::THandle<CP::TAlgorithmResult> 
    Process(const CP::TAlgorithmResult& input,
            const CP::TAlgorithmResult& input1 = CP::TAlgorithmResult::Empty,
            const CP::TAlgorithmResult& input2 = CP::TAlgorithmResult::Empty);

  CP::THitSelection ClusteredHits(CP::THitSelection& hits,bool bigest,unsigned int minPoints, double maxDist); 

  

private:

  unsigned int fminPoints;

  double fmaxDist; 
  
};



#endif
