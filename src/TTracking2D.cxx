#include "TTracking2D.hxx"
#include "HitUtilities.hxx"
#include "CompareReconObjects.hxx"
#include "CreateTrack.hxx"

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

  int CheckObjectPlane(const CP::THandle<CP::TReconBase>& obj){
    
 
    CP::THandle<CP::THitSelection> hits = obj->GetHits();
    CP::THitSelection::const_iterator beg = hits->begin();
        int plane = CP::GeomId::Captain::GetWirePlane((*beg)->GetGeomId());
    if (plane == CP::GeomId::Captain::kXPlane) {
      return 1;
    }
    if (plane == CP::GeomId::Captain::kUPlane) {
      return 2;
    }
    if (plane == CP::GeomId::Captain::kVPlane) {
      return 3;
    }

  }

std::vector<int> DeterminAxis(double LongAxis, double MajorAxis, double MinorAxis){
	    double maxAxis=std::max(LongAxis,MajorAxis);
	    maxAxis=std::max(maxAxis,MinorAxis);
	    double minAxis=std::min(LongAxis,MajorAxis);
	    minAxis=std::min(minAxis,MinorAxis);
	    int nmax=-999;
	    int nmid=-999;
	    if(maxAxis==LongAxis){
	      nmax=1;
	      if(minAxis==MajorAxis){nmid=3;}
	    else if(minAxis==MinorAxis){nmid=2;}
	    }
	    else if(maxAxis==MajorAxis){
	      nmax=2;
	      if(minAxis==MinorAxis){nmid=1;}
	      else if(minAxis==LongAxis){nmid=3;}
	    }
	    else if(maxAxis==MinorAxis){
	      nmax=3;
	      if(minAxis==MajorAxis){nmid=1;}
	      else if(minAxis==LongAxis){nmid=2;}
	    }
	    std::vector<int> ax;
	    ax.push_back(nmax);
	    ax.push_back(nmid);
	    return ax;
}

double GetAxisRation(std::vector<int> useAxis,double LongAxis, double MajorAxis, double MinorAxis){
  double ratio = -1;
  	    if(useAxis[0]==1 && useAxis[1]==2){ratio=MajorAxis/LongAxis;}
	    else if(useAxis[0]==1 && useAxis[1]==3){ratio=MinorAxis/LongAxis;}
	    else if(useAxis[0]==2 && useAxis[1]==3){ratio=MinorAxis/MajorAxis;}
	    else if(useAxis[0]==2 && useAxis[1]==1){ratio=LongAxis/MajorAxis;}
	    else if(useAxis[0]==3 && useAxis[1]==1){ratio=LongAxis/MinorAxis;}
	    else if(useAxis[0]==3 && useAxis[1]==2){ratio=MajorAxis/MinorAxis;}

	    return ratio;
}
  

CP::TTracking2D::TTracking2D()
    : TAlgorithm("TTracking2D", 
                 "Break up objects into separate hits") {

  fXZRatio=0.4;
 
}

CP::TTracking2D::~TTracking2D() { }



CP::THandle<CP::TAlgorithmResult>
CP::TTracking2D::Process(const CP::TAlgorithmResult& input,
                               const CP::TAlgorithmResult& input1,
                               const CP::TAlgorithmResult&) {

  //ReconObjects 
  CP::THandle<CP::TReconObjectContainer> inputObjects 
      = input.GetResultsContainer();
  
  //1
    CaptLog("TTracking2D Process " << GetEvent().GetContext());

    if (!inputObjects) {
        CaptError("No input objects");
        return CP::THandle<CP::TAlgorithmResult>();
    }
    
    // Create the output containers.
    CP::THandle<CP::TAlgorithmResult> result = CreateResult();
    std::unique_ptr<CP::TReconObjectContainer> 
        final(new CP::TReconObjectContainer("final"));
        std::unique_ptr<CP::TReconObjectContainer> 
        tracksX(new CP::TReconObjectContainer("tracksX"));
	    std::unique_ptr<CP::TReconObjectContainer> 
        tracksU(new CP::TReconObjectContainer("tracksU"));
	    std::unique_ptr<CP::TReconObjectContainer> 
        tracksV(new CP::TReconObjectContainer("tracksV"));

    //Copy inputObjects in "final" container
    CP::TReconObjectContainer clustersX;
    CP::TReconObjectContainer clustersU;
    CP::TReconObjectContainer clustersV;
    
    for (CP::TReconObjectContainer::iterator cl = inputObjects->begin();
         cl != inputObjects->end(); ++cl) {
        CP::THandle<CP::TReconCluster> cluster = *cl;
        if (!cluster) {
            final->push_back(*cl);
            continue;
        }
	int plane = CheckObjectPlane(cluster);
	if(plane==1){
	  std::cout<<"SEPARATING X"<<std::endl;
	  clustersX.push_back(*cl);}
	else if(plane==2){
	  clustersU.push_back(*cl);}
	else if(plane==3){
	  clustersV.push_back(*cl);}
	else{std::cout<<"PLANEDEFININGFORCLUSTERSDOESNOTWORK"<<std::endl;}
    }

    if(clustersX.size()>0)
      {
	for(CP::TReconObjectContainer::iterator cl = clustersX.begin(); cl!= clustersX.end();++cl)
	  { CP::THandle<CP::TReconCluster> cluster = *cl;
	    double LongAxis = cluster->GetLongExtent();
	    double MajorAxis = cluster->GetMajorExtent();
	    double MinorAxis = cluster->GetMinorExtent();
	    double ratio=-1;
	    std::cout<<"LongX="<<LongAxis<<"; MajorX="<<MajorAxis<<"; MinorX="<<MinorAxis<<std::endl;
	    std::vector<int> useAxis = DeterminAxis(LongAxis,MajorAxis,MinorAxis);
	    
	    ratio = GetAxisRation(useAxis,LongAxis,MajorAxis,MinorAxis);
	    std::cout<<"RatioX="<<ratio<<std::endl;
	    if(ratio<fXZRatio){
	      CP::THandle<CP::THitSelection> hft = cluster->GetHits();
	      if((*hft).size()>3){
		CP::THandle<CP::TReconTrack> track ;
		if(useAxis[0]==1){
	      track = CreateTrackFromHits("TTracking2D",
					  (*hft).begin(), (*hft).end(),cluster->GetLongAxis());
		}
		if(useAxis[0]==2){
	      track = CreateTrackFromHits("TTracking2D",
					  (*hft).begin(), (*hft).end(),cluster->GetMajorAxis());
		}
		if(useAxis[0]==3){
	      track = CreateTrackFromHits("TTracking2D",
					  (*hft).begin(), (*hft).end(),cluster->GetMinorAxis());
		}
		if(track){
                final->push_back(track);
	      tracksX->push_back(track);
		}
	      }
	    }
	  }
      }

    CaptNamedLog("TTrackingXClusters",
                 "With " << (*tracksX).size() << " tracks_xh"
                 << " from " << clustersX.size() << " x clusters");

     if(clustersU.size()>0)
      {
	for(CP::TReconObjectContainer::iterator cl = clustersU.begin(); cl!= clustersU.end();++cl)
	  { CP::THandle<CP::TReconCluster> cluster = *cl;
	    double LongAxis = cluster->GetLongExtent();
	    double MajorAxis = cluster->GetMajorExtent();
	    double MinorAxis = cluster->GetMinorExtent();
	    double ratio=-1;
	    std::cout<<"LongU="<<LongAxis<<"; MajorU="<<MajorAxis<<"; MinorU="<<MinorAxis<<std::endl;
  std::vector<int> useAxis = DeterminAxis(LongAxis,MajorAxis,MinorAxis);
	    
	    ratio = GetAxisRation(useAxis,LongAxis,MajorAxis,MinorAxis);

	    std::cout<<"RatioU="<<ratio<<std::endl;
	    if(ratio<fXZRatio){
	      CP::THandle<CP::THitSelection> hft = cluster->GetHits();
	      if((*hft).size()>3){
		CP::THandle<CP::TReconTrack> track ;
		if(useAxis[0]==1){
	      track = CreateTrackFromHits("TTracking2D",
					  (*hft).begin(), (*hft).end(),cluster->GetLongAxis());
		}
		if(useAxis[0]==2){
	      track = CreateTrackFromHits("TTracking2D",
					  (*hft).begin(), (*hft).end(),cluster->GetMajorAxis());
		}
		if(useAxis[0]==3){
	      track = CreateTrackFromHits("TTracking2D",
					  (*hft).begin(), (*hft).end(),cluster->GetMinorAxis());
		}
		if(track){
                final->push_back(track);
	      tracksU->push_back(track);
		}
	      }
	    }
	  }
      }

     CaptNamedLog("TTrackingUClusters",
                 "With " << (*tracksU).size() << " tracks_uh"
                 << " from " << clustersU.size() << " u clusters");
    if(clustersV.size()>0)
      {
	for(CP::TReconObjectContainer::iterator cl = clustersV.begin(); cl!= clustersV.end();++cl)
	  { CP::THandle<CP::TReconCluster> cluster = *cl;
	    double LongAxis = cluster->GetLongExtent();
	    double MajorAxis = cluster->GetMajorExtent();
	    double MinorAxis = cluster->GetMinorExtent();
	    double ratio=-1;
	    std::cout<<"LongV="<<LongAxis<<"; MajorV="<<MajorAxis<<"; MinorV="<<MinorAxis<<std::endl;
  std::vector<int> useAxis = DeterminAxis(LongAxis,MajorAxis,MinorAxis);
	    
	    ratio = GetAxisRation(useAxis,LongAxis,MajorAxis,MinorAxis);

	    std::cout<<"RatioV="<<ratio<<std::endl;
	    if(ratio<fXZRatio){
	      CP::THandle<CP::THitSelection> hft = cluster->GetHits();
	      if((*hft).size()>3){
		CP::THandle<CP::TReconTrack> track ;
		if(useAxis[0]==1){
	      track = CreateTrackFromHits("TTracking2D",
					  (*hft).begin(), (*hft).end(),cluster->GetLongAxis());
		}
		if(useAxis[0]==2){
	      track = CreateTrackFromHits("TTracking2D",
					  (*hft).begin(), (*hft).end(),cluster->GetMajorAxis());
		}
		if(useAxis[0]==3){
	      track = CreateTrackFromHits("TTracking2D",
					  (*hft).begin(), (*hft).end(),cluster->GetMinorAxis());
		}
		if(track){
                final->push_back(track);
	      tracksV->push_back(track);
		}
	      }
	    }
	  }
	  }
    CaptNamedLog("TTrackingVClusters",
                 "With " << (*tracksV).size() << " tracks_vh"
                 << " from " << clustersV.size() << " v clusters");
    result->AddResultsContainer(tracksX.release());
    result->AddResultsContainer(tracksU.release());
    result->AddResultsContainer(tracksV.release());
    result->AddResultsContainer(final.release());

    return result;
}
