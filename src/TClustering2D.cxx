#include "TClustering2D.hxx"
#include "HitUtilities.hxx"
#include "TPositionDensityCluster.hxx"
#include "THoughTransform.hxx"
#include "CreateCluster.hxx"
#include "CompareReconObjects.hxx"

#include <THandle.hxx>
#include <TReconTrack.hxx>
#include <TReconCluster.hxx>
#include <TCaptLog.hxx>
#include <HEPUnits.hxx>
#include <TUnitsTable.hxx>
#include <TRuntimeParameters.hxx>

#include <TWireClustering.hxx>

#include <memory>
#include <set>
#include <cmath>

#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TGraph.h"

//#define WIRE_CLUSTERING

double GetDist(double x0, double y0, std::pair<double,double> line)
{
  double a = line.first;
  double b = line.second;
  double x1=1.0;
  double x2=2.0;
  double y1=a*x1+b;
  double y2=a*x2+b;
  double d;
  double r = fabs((x2*1.0-x1*1.0)*(y1*1.0-y0*1.0)-(x1*1.0-x0*1.0)*(y2*1.0-y1*1.0));
  double l = sqrt((x2*1.0-x1*1.0)*(x2*1.0-x1*1.0)+(y2*1.0-y1*1.0)*(y2*1.0-y1*1.0));
  d=r/l;
  return d;
}

bool CompX(const CP::THandle<CP::THit>& lhs, const CP::THandle<CP::THit>& rhs){
  return lhs->GetPosition().X() < rhs->GetPosition().X();
}

CP::TClustering2D::TClustering2D()
  : TAlgorithm("TClustering2D", 
	       "Break up objects into separate hits") {

  fminPoints = CP::TRuntimeParameters::Get().GetParameterD("captRecon.clustering2D.minPoints");
  fmaxDist = CP::TRuntimeParameters::Get().GetParameterD("captRecon.clustering2D.maxDist");
  
}

CP::TClustering2D::~TClustering2D() { }

CP::THitSelection CP::TClustering2D::ClusteredHits(CP::THitSelection& hits,bool bigest,unsigned int minPoints, double maxDist){
  CP::THitSelection finalHits;
  typedef CP::TPositionDensityCluster< CP::THandle<THit> > ClusterAlgorithm;

   
  ClusterAlgorithm Clusters(minPoints,maxDist);
  Clusters.Cluster(hits.begin(), hits.end());
  int nClusters = Clusters.GetClusterCount();
  if(!bigest){
    for (int i=0; i<nClusters; ++i) {
      const ClusterAlgorithm::Points& points 
	= Clusters.GetCluster(i);
      std::copy(points.begin(),points.end(),back_inserter(finalHits));
    }
  }else{
    int maxi=-9999;
    int maxp=-9999;
    for (int i=0; i<nClusters; ++i) {
      const ClusterAlgorithm::Points& points 
	= Clusters.GetCluster(i);
      int size=points.size();
      if(size>maxp){maxp=points.size();maxi=i;}
    }
    const ClusterAlgorithm::Points& points 
      = Clusters.GetCluster(maxi);
    std::copy(points.begin(),points.end(),back_inserter(finalHits));
  }
  return finalHits;
}

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
      CP::THitSelection xHits_cl;
      
      xHits_cl=ClusteredHits(xHits,false,5,25);
      bool failHough=false;
      typedef CP::THoughTrans< CP::THandle<THit> > HoughAlgorithm;
      
      while(!failHough){
	HoughAlgorithm Hough(180,5000);

	std::vector<std::pair<double,double>> points;

	for(CP::THitSelection::iterator h=xHits_cl.begin();h!=xHits_cl.end();++h){
	  double x = (*h)->GetPosition().X();
	  points.push_back(std::make_pair(x,(*h)->GetPosition().Z()));
	}
	
	Hough.HoughTransform(points.begin(),points.end());

	
	std::pair<double,double> lineParam=Hough.GetLineParam(10);

	if(lineParam.first==-9999 || lineParam.second==-9999){
	  failHough=true;
	  continue;
	}
	
	CP::THitSelection LineHits;
	CP::THitSelection FinalLine;

	for(CP::THitSelection::iterator h=xHits_cl.begin();h!=xHits_cl.end();++h){
	  double dist = GetDist((*h)->GetPosition().X(),(*h)->GetPosition().Z(),lineParam);

	  if(dist<15*unit::mm) LineHits.push_back(*h);
	}

	if(LineHits.size()>2){
	  std::sort(LineHits.begin(),LineHits.end(),CompX);
	  FinalLine=CP::TClustering2D::ClusteredHits(LineHits,true,3,60);
	  //	std::copy(LineHits.begin(),LineHits.end(),back_inserter(FinalLine));

	  CP::THandle<CP::TReconCluster> cluster = CreateCluster("cluster_xh",FinalLine.begin(),FinalLine.end());
	  if(cluster){
#ifndef WIRE_CLUSTERING
	    xclusters->push_back(cluster);
	    final->push_back(cluster);
#endif
	  }
	}
	
	if(FinalLine.size()>0){
	  for(CP::THitSelection::iterator h=FinalLine.begin();h!=FinalLine.end();++h){
	    xHits_cl.erase(std::remove(xHits_cl.begin(),xHits_cl.end(),*h),xHits_cl.end());
	  }
	}else{failHough=true;}
	if(xHits_cl.size()<5)
	  failHough=true;
	points.clear();
	  
      }

#ifdef WIRE_CLUSTERING
      CP::TWireClustering wireClusters;
      for (auto cluster:wireClusters.WireCluster(0,xHits,10,6,50)) {
	    xclusters->push_back(cluster);
	    final->push_back(cluster);
      }
#endif
      std::cout<<"NClusters="<<xclusters->size()<<std::endl;
     
      CaptNamedLog("TClusterXHits",
		   "With " << xclusters->size() << " clusters_xh");
    
    }

  if(uHits.size()>0)
    {

      CP::THitSelection uHits_cl;

      uHits_cl=ClusteredHits(uHits,false,3,25);
      //		std::copy(uHits.begin(),uHits.end(),back_inserter(uHits_cl));
      bool failHough=false;
      typedef CP::THoughTrans< CP::THandle<THit> > HoughAlgorithm;

      while(!failHough){
	HoughAlgorithm Hough(180,7000);

	
	std::vector<std::pair<double,double>> points;

	for(CP::THitSelection::iterator h=uHits_cl.begin();h!=uHits_cl.end();++h){
	  double sign = 1.0;
	  if((*h)->GetPosition().X()<0)sign=-1.0;
	  double x=sign*sqrt((*h)->GetPosition().X()*(*h)->GetPosition().X()+(*h)->GetPosition().Y()*(*h)->GetPosition().Y());//(sin(30*3.14159265/180.0)*sin(30*3.14159265/180.0));
	  points.push_back(std::make_pair(x,(*h)->GetPosition().Z()));
	}
	
	Hough.HoughTransform(points.begin(),points.end());

	
	std::pair<double,double> lineParam=Hough.GetLineParam(20);

	if(lineParam.first==-9999 || lineParam.second==-9999){
	  failHough=true;
	  continue;
	}
	
	CP::THitSelection LineHits;
	CP::THitSelection FinalLine;
	for(CP::THitSelection::iterator h=uHits_cl.begin();h!=uHits_cl.end();++h){
	  double sign = 1.0;
	  if((*h)->GetPosition().X()<0)sign=-1.0;
	  double x=sign*sqrt((*h)->GetPosition().X()*(*h)->GetPosition().X()+(*h)->GetPosition().Y()*(*h)->GetPosition().Y());
	  double dist = GetDist(x,(*h)->GetPosition().Z(),lineParam);

	  if(dist<20) LineHits.push_back(*h);
	}

	if(LineHits.size()>2){
	  std::sort(LineHits.begin(),LineHits.end(),CompX);
	  FinalLine=CP::TClustering2D::ClusteredHits(LineHits,true,2,60);
	  //	std::copy(LineHits.begin(),LineHits.end(),back_inserter(FinalLine));

	  CP::THandle<CP::TReconCluster> cluster = CreateCluster("cluster_uh",FinalLine.begin(),FinalLine.end());
	  if(cluster){
#ifndef WIRE_CLUSTERING
	    uclusters->push_back(cluster);
	    final->push_back(cluster);
#endif
	  }
	}
	
	if(FinalLine.size()>0){
	  for(CP::THitSelection::iterator h=FinalLine.begin();h!=FinalLine.end();++h){
	    uHits_cl.erase(std::remove(uHits_cl.begin(),uHits_cl.end(),*h),uHits_cl.end());
	  }
	}else{failHough=true;}
	if(uHits_cl.size()<5)
	  failHough=true;
	points.clear();
	  
      }

#ifdef WIRE_CLUSTERING	
      CP::TWireClustering wireClusters;
      for (auto cluster:wireClusters.WireCluster(2,uHits,10,6,50)) {
	uclusters->push_back(cluster);
	final->push_back(cluster);
      }
#endif
      std::cout<<"NClusters="<<uclusters->size()<<std::endl;
     
      CaptNamedLog("TClusterUHits",
		   "With " << uclusters->size() << " clusters_uh");
    
    }

     
  if(vHits.size()>0)
    {
      CP::THitSelection vHits_cl;

      vHits_cl=ClusteredHits(vHits,false,5,25);
      //		std::copy(vHits.begin(),vHits.end(),back_inserter(vHits_cl));
      bool failHough=false;
      typedef CP::THoughTrans< CP::THandle<THit> > HoughAlgorithm;

      while(!failHough){
	HoughAlgorithm Hough(180,3000);

	
	std::vector<std::pair<double,double>> points;

	for(CP::THitSelection::iterator h=vHits_cl.begin();h!=vHits_cl.end();++h){
	  double sign = 1.0;
	  if((*h)->GetPosition().X()<0)sign=-1.0;
	  double x=sign*sqrt((*h)->GetPosition().X()*(*h)->GetPosition().X()+(*h)->GetPosition().Y()*(*h)->GetPosition().Y());//(sin(30*3.14159265/180.0)*sin(30*3.14159265/180.0));
	  points.push_back(std::make_pair(x,(*h)->GetPosition().Z()));
	}
	
	Hough.HoughTransform(points.begin(),points.end());

	
	std::pair<double,double> lineParam=Hough.GetLineParam(20);

	if(lineParam.first==-9999 || lineParam.second==-9999){
	  failHough=true;
	  continue;
	}
	
	CP::THitSelection LineHits;
	for(CP::THitSelection::iterator h=vHits_cl.begin();h!=vHits_cl.end();++h){
	  double sign = 1.0;
	  if((*h)->GetPosition().X()<0)sign=-1.0;
	  double x=sign*sqrt((*h)->GetPosition().X()*(*h)->GetPosition().X()+(*h)->GetPosition().Y()*(*h)->GetPosition().Y());
	  double dist = GetDist(x,(*h)->GetPosition().Z(),lineParam);

	  if(dist<10) LineHits.push_back(*h);
	}

	CP::THitSelection FinalLine;
	
	if(LineHits.size()>2){
	  std::sort(LineHits.begin(),LineHits.end(),CompX);

	  FinalLine=CP::TClustering2D::ClusteredHits(LineHits,true,2,60);
	  //	std::copy(LineHits.begin(),LineHits.end(),back_inserter(FinalLine));

	  CP::THandle<CP::TReconCluster> cluster = CreateCluster("cluster_vh",FinalLine.begin(),FinalLine.end());
	  if(cluster){
#ifndef WIRE_CLUSTERING
	    vclusters->push_back(cluster);
	    final->push_back(cluster);
#endif
	  }
	}
	
	if(FinalLine.size()>0){
	  for(CP::THitSelection::iterator h=FinalLine.begin();h!=FinalLine.end();++h){
	    vHits_cl.erase(std::remove(vHits_cl.begin(),vHits_cl.end(),*h),vHits_cl.end());
	  }
	}else{failHough=true;}
	if(vHits_cl.size()<5)
	  failHough=true;
	points.clear();
	
      }

#ifdef WIRE_CLUSTERING
      CP::TWireClustering wireClusters;
      for (auto cluster:wireClusters.WireCluster(1,vHits,10,6,50)) {
	vclusters->push_back(cluster);
	final->push_back(cluster);
      }
#endif
      std::cout<<"NClusters="<<vclusters->size()<<std::endl;
	
     
      CaptNamedLog("TClusterVHits",
		   "With " << vclusters->size() << " clusters_vh");
    
    }
    
  std::unique_ptr<TH2F> HitsX(new TH2F("HitsForX","HitsForX",340,0,340,9600,0,9600));
  std::unique_ptr<TH2F> HitsU(new TH2F("HitsForU","HitsForU",340,0,340,9600,0,9600));
  std::unique_ptr<TH2F>  HitsV(new TH2F("HitsForV","HitsForV",340,0,340,9600,0,9600));
  
  if(xclusters->size()>0){
    int coll=1;
    for(CP::TReconObjectContainer::iterator it = xclusters->begin();it!=xclusters->end();++it){
      CP::THandle<CP::THitSelection> hits = (*it)->GetHits();
      if(hits){
	for(CP::THitSelection::iterator h = hits->begin();h!=hits->end();++h){
	  double ht=((*h)->GetConstituent()->GetTime()+1.6*unit::ms)/(500*unit::ns);
	  int j = HitsX->GetYaxis()->FindFixBin(ht);
	  double hw=CP::GeomId::Captain::GetWireNumber((*h)->GetGeomId());
	  int i = HitsX->GetXaxis()->FindFixBin(hw);
	  HitsX->SetBinContent(i,j,coll*10);
	}

      }
      ++coll;
    }

    HitsX->Draw("COLZ");
    gPad->Print("plots/XHits_cl.C");
  }

  if(uclusters->size()>0){
    int coll=1;
    for(CP::TReconObjectContainer::iterator it = uclusters->begin();it!=uclusters->end();++it){
      CP::THandle<CP::THitSelection> hits = (*it)->GetHits();
      if(hits){
	for(CP::THitSelection::iterator h = hits->begin();h!=hits->end();++h){
	  double ht=((*h)->GetConstituent()->GetTime()+1.6*unit::ms)/(500*unit::ns);
	  int j = HitsU->GetYaxis()->FindFixBin(ht);
	  double hw=CP::GeomId::Captain::GetWireNumber((*h)->GetGeomId());
	  int i = HitsU->GetXaxis()->FindFixBin(hw);
	  HitsU->SetBinContent(i,j,coll*10);
	}
      }
      ++coll;
    }
    
    HitsU->Draw("COLZ");
    gPad->Print("plots/UHits_cl.C");
  }

  if(vclusters->size()>0){
    int coll=1;
    for(CP::TReconObjectContainer::iterator it = vclusters->begin();it!=vclusters->end();++it){
      CP::THandle<CP::THitSelection> hits = (*it)->GetHits();
      if(hits){
	for(CP::THitSelection::iterator h = hits->begin();h!=hits->end();++h){
	  double ht=((*h)->GetConstituent()->GetTime()+1.6*unit::ms)/(500*unit::ns);
	  int j = HitsV->GetYaxis()->FindFixBin(ht);
	  double hw=CP::GeomId::Captain::GetWireNumber((*h)->GetGeomId());
	  int i = HitsV->GetXaxis()->FindFixBin(hw);
	  HitsV->SetBinContent(i,j,coll*10);
	}
      }
      coll++;
    }
 
    HitsV->Draw("COLZ");
    gPad->Print("plots/VHits_cl.C");
  }    
   

  result->AddResultsContainer(xclusters.release());
  result->AddResultsContainer(uclusters.release());
  result->AddResultsContainer(vclusters.release());
  result->AddResultsContainer(final.release());

  return result;
}
