#include "TClustering2D.hxx"
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

#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"

CP::TClustering2D::TClustering2D()
    : TAlgorithm("TClustering2D", 
                 "Break up objects into separate hits") {

  fminPoints = CP::TRuntimeParameters::Get().GetParameterD("captRecon.clustering2D.minPoints");
  fmaxDist = CP::TRuntimeParameters::Get().GetParameterD("captRecon.clustering2D.maxDist");
  
 
}

CP::TClustering2D::~TClustering2D() { }



CP::THandle<CP::TAlgorithmResult>
CP::TClustering2D::Process(const CP::TAlgorithmResult& input,
                               const CP::TAlgorithmResult&,
                               const CP::TAlgorithmResult&) {

  //ReconObjects 
  CP::THandle<CP::TReconObjectContainer> inputObjects 
      = input.GetResultsContainer();
  
  //all Hits in the event
  CP::THandle<CP::THitSelection> inputObjects1 
    = input.GetHits();
  //1
    CaptLog("TClustering2D Process " << GetEvent().GetContext());

    if (!inputObjects || !inputObjects1) {
        CaptError("No input objects");
        return CP::THandle<CP::TAlgorithmResult>();
    }
    CP::THitSelection xHits;
    CP::THitSelection vHits;
    CP::THitSelection uHits;
    for (CP::THitSelection::iterator h = inputObjects1->begin(); 
         h != inputObjects1->end(); ++h) {

        int plane = CP::GeomId::Captain::GetWirePlane((*h)->GetGeomId());

        if (plane == CP::GeomId::Captain::kXPlane) {

	   xHits.push_back(*h);
	   // std::cout<<"xHit"<<"; X="<<(*h)->GetPosition().X()<<"; Y="<<(*h)->GetPosition().Y()<<"; Z="<<(*h)->GetPosition().Z()<<"; Unc.Z="<<(*h)->GetUncertainty().Z()<<std::endl;
        }
        else if (plane == CP::GeomId::Captain::kVPlane) {
	  vHits.push_back(*h);
	  //  std::cout<<"uHit"<<"; X="<<(*h)->GetPosition().X()<<"; Y="<<(*h)->GetPosition().Y()<<"; Z="<<(*h)->GetPosition().Z()<<"; Unc.Z="<<(*h)->GetUncertainty().Z()<<std::endl;
        }
        else if (plane == CP::GeomId::Captain::kUPlane) {
	  uHits.push_back(*h);
	  // std::cout<<"vHit"<<"; X="<<(*h)->GetPosition().X()<<"; Y="<<(*h)->GetPosition().Y()<<"; Z="<<(*h)->GetPosition().Z()<<"; Unc.Z="<<(*h)->GetUncertainty().Z()<<std::endl;
        }
        else {
            CaptError("Invalid wire plane");
        }
    }
 
    // Create the output containers.
    CP::THandle<CP::TAlgorithmResult> result = CreateResult();
    std::unique_ptr<CP::TReconObjectContainer> 
        final(new CP::TReconObjectContainer("final"));
     std::unique_ptr<CP::TReconObjectContainer> xclusters(
        new CP::TReconObjectContainer("xclusters"));
      std::unique_ptr<CP::TReconObjectContainer> uclusters(
        new CP::TReconObjectContainer("uclusters"));
       std::unique_ptr<CP::TReconObjectContainer> vclusters(
        new CP::TReconObjectContainer("vclusters"));

    //Copy inputObjects in "final" container
    for (CP::TReconObjectContainer::iterator it = inputObjects->begin();
         it != inputObjects->end(); ++it)
      {
	final->push_back(*it);
      }

    //If unused hits exist, cluster them by position and save in "final" container    
    if(xHits.size()>0)
      {
    typedef CP::TPositionDensityCluster< CP::THandle<THit> > ClusterAlgorithm;

    // Find the big clusters.
    ClusterAlgorithm Clusters(fminPoints,fmaxDist);
    Clusters.Cluster(xHits.begin(), xHits.end());
    
    int nClusters = Clusters.GetClusterCount();
    CaptNamedLog("TClusterXHits",
                 "With " << nClusters << " clusters_xh"
                 << " from " << xHits.size() << " x hits");
    for (int i=0; i<nClusters; ++i) {
        const ClusterAlgorithm::Points& points 
            = Clusters.GetCluster(i);
        CP::THandle<CP::TReconCluster> cluster
            = CreateCluster("cluster_xh",points.begin(),points.end());
        final->push_back(cluster);
	xclusters->push_back(cluster);
    }
   std::unique_ptr<CP::THitSelection> xused(new CP::THitSelection("xused"));
   std::unique_ptr<CP::THitSelection> xunused(new CP::THitSelection("xunused"));

   std::copy(xHits.begin(), xHits.end(), std::back_inserter(*xused));
   xunused->reserve(xHits.size());
        std::copy(xHits.begin(), xHits.end(), 
                  std::back_inserter(*xunused));
        CP::hits::Subtract(*xunused,*xused);

	 result->AddHits(xunused.release());
        result->AddHits(xused.release());
      }

     if(uHits.size()>0)
      {
    typedef CP::TPositionDensityCluster< CP::THandle<THit> > ClusterAlgorithm;

    // Find the big clusters.
    ClusterAlgorithm Clusters(fminPoints,fmaxDist);
    Clusters.Cluster(uHits.begin(), uHits.end());
    //for(std::size_t i=0;i<uHits.size();++i)
    // std::cout<<(*uHits[i]).GetPosition().Z()<<std::endl;
    int nClusters = Clusters.GetClusterCount();
    CaptNamedLog("TClusterUHits",
                 "With " << nClusters << " clusters_uh"
                 << " from " << uHits.size() << " u hits");
    for (int i=0; i<nClusters; ++i) {
        const ClusterAlgorithm::Points& points 
            = Clusters.GetCluster(i);
        CP::THandle<CP::TReconCluster> cluster
            = CreateCluster("cluster_uh",points.begin(),points.end());
        final->push_back(cluster);
	uclusters->push_back(cluster);
    }
    
   std::unique_ptr<CP::THitSelection> uused(new CP::THitSelection("uused"));
   std::unique_ptr<CP::THitSelection> uunused(new CP::THitSelection("uunused"));

   std::copy(uHits.begin(), uHits.end(), std::back_inserter(*uused));
   uunused->reserve(uHits.size());
        std::copy(uHits.begin(), uHits.end(), 
                  std::back_inserter(*uunused));
        CP::hits::Subtract(*uunused,*uused);

	 result->AddHits(uunused.release());
        result->AddHits(uused.release());
      }

     
     if(vHits.size()>0)
      {
    typedef CP::TPositionDensityCluster< CP::THandle<THit> > ClusterAlgorithm;

    // Find the big clusters.
    ClusterAlgorithm Clusters(fminPoints,fmaxDist);
    Clusters.Cluster(vHits.begin(), vHits.end());
    
    int nClusters = Clusters.GetClusterCount();
    CaptNamedLog("TClusterVHits",
                 "With " << nClusters << " clusters_vh"
                 << " from " << vHits.size() << " v hits");
    for (int i=0; i<nClusters; ++i) {
        const ClusterAlgorithm::Points& points 
            = Clusters.GetCluster(i);
        CP::THandle<CP::TReconCluster> cluster
            = CreateCluster("cluster_vh",points.begin(),points.end());
        final->push_back(cluster);
	vclusters->push_back(cluster);
    }
    
   std::unique_ptr<CP::THitSelection> vused(new CP::THitSelection("vused"));
   std::unique_ptr<CP::THitSelection> vunused(new CP::THitSelection("vunused"));

   std::copy(vHits.begin(), vHits.end(), std::back_inserter(*vused));
   vunused->reserve(vHits.size());
        std::copy(vHits.begin(), vHits.end(), 
                  std::back_inserter(*vunused));
        CP::hits::Subtract(*vunused,*vused);

	 result->AddHits(vunused.release());
        result->AddHits(vused.release());
      } 
    
    TH2F* HitsX = new TH2F("HitsForX","HitsForX",340,0,340,9600,0,9600);
    TH2F* HitsU = new TH2F("HitsForU","HitsForU",340,0,340,9600,0,9600);
    TH2F* HitsV = new TH2F("HitsForV","HitsForV",340,0,340,9600,0,9600);
    if(xclusters->size()>0){
      for(CP::TReconObjectContainer::iterator it = xclusters->begin();it!=xclusters->end();++it){
	CP::THandle<CP::THitSelection> hits = (*it)->GetHits();
	if(hits){
	  for(CP::THitSelection::iterator h = hits->begin();h!=hits->end();++h){
	    double ht=((*h)->GetConstituent()->GetTime()+1.6*unit::ms)/(500*unit::ns);
	    double hw=CP::GeomId::Captain::GetWireNumber((*h)->GetGeomId());
	    HitsX->Fill(hw,ht);
	  }
	}
      }

      HitsX->Draw();
      gPad->Print("plots/XHits_cl.C");
    }

      if(uclusters->size()>0){
      for(CP::TReconObjectContainer::iterator it = uclusters->begin();it!=uclusters->end();++it){
	CP::THandle<CP::THitSelection> hits = (*it)->GetHits();
	if(hits){
	  for(CP::THitSelection::iterator h = hits->begin();h!=hits->end();++h){
	    double ht=((*h)->GetConstituent()->GetTime()+1.6*unit::ms)/(500*unit::ns);
	    double hw=CP::GeomId::Captain::GetWireNumber((*h)->GetGeomId());
	    HitsU->Fill(hw,ht);

	  }
	}
      }

      HitsU->Draw();
      gPad->Print("plots/UHits_cl.C");
    }

        if(vclusters->size()>0){
      for(CP::TReconObjectContainer::iterator it = vclusters->begin();it!=vclusters->end();++it){
	CP::THandle<CP::THitSelection> hits = (*it)->GetHits();
	if(hits){
	  for(CP::THitSelection::iterator h = hits->begin();h!=hits->end();++h){
	    double ht=((*h)->GetConstituent()->GetTime()+1.6*unit::ms)/(500*unit::ns);
	    double hw=CP::GeomId::Captain::GetWireNumber((*h)->GetGeomId());
	    HitsV->Fill(hw,ht);

	  }
	}
      }
 
      HitsV->Draw();
      gPad->Print("plots/VHits_cl.C");
    }

     
   

    result->AddResultsContainer(xclusters.release());
     result->AddResultsContainer(uclusters.release());
    result->AddResultsContainer(vclusters.release());
    result->AddResultsContainer(final.release());

    return result;
}
