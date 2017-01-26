#include "TClusterUnusedHits.hxx"
#include "HitUtilities.hxx"
#include "TPositionDensityCluster.hxx"
#include "CreateCluster.hxx"
#include "CompareReconObjects.hxx"

#include <THandle.hxx>
#include <TReconTrack.hxx>
#include <TReconCluster.hxx>
#include <TCaptLog.hxx>
#include <HEPUnits.hxx>
#include <TUnitsTable.hxx>
#include <TRuntimeParameters.hxx>
 
#include <memory>
#include <set>
#include <cmath>

CP::TClusterUnusedHits::TClusterUnusedHits()
    : TAlgorithm("TClusterUnusedHits", 
                 "Break up objects into separate hits") {

  fminPoints = CP::TRuntimeParameters::Get().GetParameterD("captRecon.clusterUnusedHits.minPoints");
  fmaxDist = CP::TRuntimeParameters::Get().GetParameterD("captRecon.clusterUnusedHits.maxDist");
 
}

CP::TClusterUnusedHits::~TClusterUnusedHits() { }

CP::THandle<CP::TAlgorithmResult>
CP::TClusterUnusedHits::Process(const CP::TAlgorithmResult& input,
                               const CP::TAlgorithmResult& input1,
                               const CP::TAlgorithmResult&) {

  //ReconObjects 
  CP::THandle<CP::TReconObjectContainer> inputObjects 
      = input.GetResultsContainer();
  
  //all Hits in the event
  CP::THandle<CP::THitSelection> inputObjects1 
    = input1.GetHits();
  //1
    CaptLog("TClusterUnusedHits Process " << GetEvent().GetContext());

    if (!inputObjects || !inputObjects1) {
        CaptError("No input objects");
        return CP::THandle<CP::TAlgorithmResult>();
    }

    //Find Hits that were not used by ReconObjects in inputObjects.
    
     std::unique_ptr<CP::THitSelection> 
            used(new CP::THitSelection("used_uh"));
        std::unique_ptr<CP::THitSelection> 
            unused(new CP::THitSelection("unused_uh"));
    
     CP::THandle<CP::THitSelection> hits 
            = CP::hits::ReconHits(inputObjects->begin(), inputObjects->end());
        if (hits) {
            std::copy(hits->begin(), hits->end(), std::back_inserter(*used));
        }
        
        unused->reserve(inputObjects1->size());
        std::copy(inputObjects1->begin(), inputObjects1->end(), 
                  std::back_inserter(*unused));
        CP::hits::Subtract(*unused,*used);


    // Create the output containers.
    CP::THandle<CP::TAlgorithmResult> result = CreateResult();
    std::unique_ptr<CP::TReconObjectContainer> 
        final(new CP::TReconObjectContainer("final"));

    //Copy inputObjects in "final" container
    for (CP::TReconObjectContainer::iterator it = inputObjects->begin();
         it != inputObjects->end(); ++it)
      {
	final->push_back(*it);
      }

    //If unused hits exist, cluster them by position and save in "final" container
    
    if(unused)
      {
    typedef CP::TPositionDensityCluster< CP::THandle<THit> > ClusterAlgorithm;

    // Find the big clusters.
    ClusterAlgorithm Clusters(fminPoints,fmaxDist);
    Clusters.Cluster(unused->begin(), unused->end());
    
    int nClusters = Clusters.GetClusterCount();
    CaptNamedLog("TClusterUnusedHits",
                 "With " << nClusters << " clusters_uh"
                 << " from " << unused->size() << " unused hits");
    for (int i=0; i<nClusters; ++i) {
        const ClusterAlgorithm::Points& points 
            = Clusters.GetCluster(i);
        CP::THandle<CP::TReconCluster> cluster
            = CreateCluster("cluster_uh",points.begin(),points.end());
        final->push_back(cluster);
    }
    
      } 
  
    result->AddResultsContainer(final.release());

    return result;
}
