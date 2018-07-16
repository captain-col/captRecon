#ifndef TTracking3D_hxx_seen
#define TTracking3D_hxx_seen

#include <TAlgorithm.hxx>
#include <TAlgorithmResult.hxx>
#include <THandle.hxx>
#include <TReconTrack.hxx>
#include <TReconCluster.hxx>

#include "TH2F.h"
#include "TGraph.h"
#include "TMultiGraph.h"

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

  bool Assemble3DTrack( CP::THandle<CP::TReconCluster> trackX, CP::THandle<CP::TReconCluster> trackU, CP::THandle<CP::TReconCluster> trackV,CP::TReconObjectContainer& match3,int trackNum);

  void FindTrackCandidates(CP::TReconObjectContainer& tracksX,CP::TReconObjectContainer& tracksU,CP::TReconObjectContainer& tracksV,CP::TReconObjectContainer& match3,CP::TReconObjectContainer& match2);

  bool Assemble2DTrack( CP::THandle<CP::TReconCluster> trackX, CP::THandle<CP::TReconCluster> trackU,CP::TReconObjectContainer& match2,int trackNum,int planeComb);
  

private:

  TMultiGraph* fHitsX;// = new TH2F("HitsForX","HitsForX",340,0,340,9600,0,9600);
  TMultiGraph* fHitsU;// = new TH2F("HitsForU","HitsForU",340,0,340,9600,0,9600));
  TMultiGraph* fHitsV;//(new TH2F("HitsForV","HitsForV",340,0,340,9600,0,9600));


};



#endif
