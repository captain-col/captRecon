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
#include "TF1.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TGraph.h"
#include "TVectorD.h"

//#define WIRE_CLUSTERING

double GetDist(double x0, double y0, std::pair<double,double> line)
{
  double a = line.first;
  double b = line.second;
  double d;
  d=fabs(a*x0-y0+b)/sqrt(a*a+1);

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

for(CP::THitSelection::iterator h=hits.begin();h!=hits.end();++h){
	 
	  double ht=((*h)->GetTime()+1.6*unit::ms)/(500*unit::ns);
	  double hw=CP::GeomId::Captain::GetWireNumber((*h)->GetGeomId());
       
	  //  std::cout<<"t="<<ht<<"; w="<<hw<<std::endl;

	  //  std::cout<<"X="<<(*h)->GetPosition().X()<<"; Y="<<(*h)->GetPosition().Y()<<"; Z="<<(*h)->GetPosition().Z()<<std::endl;

	  if(h<hits.end()-1){
	    //   std::cout<<"distance="<< ((*h)->GetPosition()-(*(h+1))->GetPosition()).Mag()<<std::endl;
				 }
 }

  
  ClusterAlgorithm Clusters(minPoints,maxDist);
  Clusters.Cluster(hits.begin(), hits.end());
  int nClusters = Clusters.GetClusterCount();
  // std::cout<<"nclusters="<<nClusters<<std::endl;
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
    if(nClusters>0)
    std::copy(points.begin(),points.end(),back_inserter(finalHits));
  }
  return finalHits;
}

CP::THandle<CP::TAlgorithmResult>
CP::TClustering2D::Process(const CP::TAlgorithmResult& input,
			   const CP::TAlgorithmResult&,
			   const CP::TAlgorithmResult&) {

  int finallinecut=0;

  TH1F* distance = new TH1F("dist","dist",30,0,30);
  TH2F* distance2D = new TH2F("dist2D","dis2D",30,0,30,20,0,20);
  TH1F* chi = new TH1F("chi","chi",1000,0,1000);
  TH2F* chi_angle = new TH2F("chi_angle","chi_angle",1000,0,1000,180,0,180);
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
  std::unique_ptr<CP::THitSelection> hitsAll(new CP::THitSelection(
                                                     "hitsAll"));
  std::unique_ptr<CP::TReconObjectContainer> 
    final(new CP::TReconObjectContainer("final"));
  std::unique_ptr<CP::TReconObjectContainer> xclusters(
						       new CP::TReconObjectContainer("xclusters"));
  std::unique_ptr<CP::TReconObjectContainer> uclusters(
						       new CP::TReconObjectContainer("uclusters"));
  std::unique_ptr<CP::TReconObjectContainer> vclusters(
						       new CP::TReconObjectContainer("vclusters"));

  //Copy inputObjects in "final" container
  for (CP::THitSelection::iterator it = inputObjects1->begin();
       it != inputObjects1->end(); ++it)
    {
      CP::THandle<CP::TReconHit> hit = *it;
      if(hit){
	 hitsAll->push_back(hit);
      }
    }
  
  //If unused hits exist, cluster them by position and save in "final" container

  if(xHits.size()>0)
    {
      CP::THitSelection xHits_cl;
      
      // xHits_cl=ClusteredHits(xHits,false,5,25);
      std::copy(xHits.begin(),xHits.end(),back_inserter(xHits_cl));
      bool failHough=false;
      typedef CP::THoughTrans< CP::THandle<THit> > HoughAlgorithm;
      
      while(!failHough){
	HoughAlgorithm Hough(2*180,10000);

	std::vector<std::pair<double,double>> points;

	for(CP::THitSelection::iterator h=xHits_cl.begin();h!=xHits_cl.end();++h){
	  double x = (*h)->GetPosition().X();
	  points.push_back(std::make_pair(x,(*h)->GetPosition().Z()));
	}
	
	Hough.HoughTransform(points.begin(),points.end());


	
	std::pair<double,double> lineParam=Hough.GetLineParam(3);

	
	if(lineParam.first==-9999 || lineParam.second==-9999){
	  failHough=true;
	  continue;
	}
	
	CP::THitSelection LineHits;
	CP::THitSelection FinalLine;

	for(CP::THitSelection::iterator h=xHits_cl.begin();h!=xHits_cl.end();++h){
	  double dist = GetDist((*h)->GetPosition().X(),(*h)->GetPosition().Z(),lineParam);
	  // std::cout<<"dist = "<<dist<<std::endl;
	  //	  distance->Fill(dist);
	  if(dist<4*unit::mm) LineHits.push_back(*h);
	}
	//	std::cout<<"lineHits="<<LineHits.size()<<std::endl;
	if(LineHits.size()>2){
	  std::sort(LineHits.begin(),LineHits.end(),CompX);
	   
	  FinalLine=CP::TClustering2D::ClusteredHits(LineHits,true,1,53*unit::mm);
	  // std::cout<<"finalLine="<<FinalLine.size()<<std::endl;
	  //	std::copy(LineHits.begin(),LineHits.end(),back_inserter(FinalLine));
        
	  if(FinalLine.size()>2){
	  CP::THandle<CP::TReconCluster> cluster = CreateCluster("cluster_xh",FinalLine.begin(),FinalLine.end());
	  if(cluster){
#ifndef WIRE_CLUSTERING

	     CP::THandle<CP::THitSelection> hits = cluster->GetHits();
	     	double avgDist=0;
      if(hits){

	for(CP::THitSelection::iterator h = hits->begin();h!=(hits->end()-1);++h){avgDist+= ((*h)->GetPosition()- (*(h+1))->GetPosition()).Mag();
	}
      }
	if(cluster->GetHits()->size()>2 ){
	    xclusters->push_back(cluster);
	    final->push_back(cluster);
	    }
#endif
	  }
	  }
	}
	
	if(FinalLine.size()>finallinecut){
	  for(CP::THitSelection::iterator h=FinalLine.begin();h!=FinalLine.end();++h){
	    xHits_cl.erase(std::remove(xHits_cl.begin(),xHits_cl.end(),*h),xHits_cl.end());
	  }
	}
if(FinalLine.size()==0 && LineHits.size()>0){
	    for(CP::THitSelection::iterator h=LineHits.begin();h!=LineHits.end();++h){
	    xHits_cl.erase(std::remove(xHits_cl.begin(),xHits_cl.end(),*h),xHits_cl.end());
	  }
	  }
  if(FinalLine.size()==0 && LineHits.size()<3){failHough=true;}
	
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

      if(xclusters->size()>0){
    for(CP::TReconObjectContainer::iterator it = xclusters->begin();it!=xclusters->end();++it){
      CP::THandle<CP::THitSelection> hits = (*it)->GetHits();
      if(hits){
	double avgDist=0;
	for(CP::THitSelection::iterator h = hits->begin();h!=(hits->end()-1);++h){avgDist+= ((*h)->GetPosition()- (*(h+1))->GetPosition()).Mag();
	}
	avgDist/=((double)hits->size()-1);
	distance2D->Fill(avgDist,(int)hits->size());
      }
    }
      }

      
      std::cout<<"NClusters="<<xclusters->size()<<std::endl;
     
      CaptNamedLog("TClusterXHits",
		   "With " << xclusters->size() << " clusters_xh");
    
    }

  if(uHits.size()>0)
    {

      CP::THitSelection uHits_cl;

      
      
      	std::copy(uHits.begin(),uHits.end(),back_inserter(uHits_cl));

      
	
      bool failHough=false;
      typedef CP::THoughTrans< CP::THandle<THit> > HoughAlgorithm;

      while(!failHough){
	HoughAlgorithm Hough(2*180,10000);

	
	std::vector<std::pair<double,double>> points;

	for(CP::THitSelection::iterator h=uHits_cl.begin();h!=uHits_cl.end();++h){
	  double sign = 1.0;
	  if((*h)->GetPosition().X()<0)sign=-1.0;
	  double x=sign*sqrt((*h)->GetPosition().X()*(*h)->GetPosition().X()+(*h)->GetPosition().Y()*(*h)->GetPosition().Y());//(sin(30*3.14159265/180.0)*sin(30*3.14159265/180.0));
	  points.push_back(std::make_pair(x,(*h)->GetPosition().Z()));
	}
	
	Hough.HoughTransform(points.begin(),points.end());

	
	std::pair<double,double> lineParam=Hough.GetLineParam(3);
	

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
	  double ht=((*h)->GetTime()+1.6*unit::ms)/(500*unit::ns);
	  double hw=CP::GeomId::Captain::GetWireNumber((*h)->GetGeomId());
	  //  if(hw>33 && hw<50)
	    // std::cout<<"t="<<ht<<"; w="<<hw<<"; dist="<<dist<<std::endl;
	  if(dist<4*unit::mm) LineHits.push_back(*h);
	}
	//	std::cout<<"lihesize="<<LineHits.size()<<std::endl;
	if(LineHits.size()>2){
	  std::sort(LineHits.begin(),LineHits.end(),CompX);
	  FinalLine=CP::TClustering2D::ClusteredHits(LineHits,true,1,53*unit::mm);
	  // std::cout<<"filanline="<<FinalLine.size()<<std::endl;
	  if(FinalLine.size()>finallinecut){
	  CP::THandle<CP::TReconCluster> cluster = CreateCluster("cluster_uh",FinalLine.begin(),FinalLine.end());
	  if(cluster){
#ifndef WIRE_CLUSTERING
	   CP::THandle<CP::THitSelection> hits = cluster->GetHits();
	     	double avgDist=0;
      if(hits){

	for(CP::THitSelection::iterator h = hits->begin();h!=(hits->end()-1);++h){avgDist+= ((*h)->GetPosition()- (*(h+1))->GetPosition()).Mag();
	}
      }
	if(cluster->GetHits()->size()>2 ){
	    uclusters->push_back(cluster);
	    final->push_back(cluster);
	    }
#endif
	  }
	  }
	}
	
	if(FinalLine.size()>0){
	  for(CP::THitSelection::iterator h=FinalLine.begin();h!=FinalLine.end();++h){
	    uHits_cl.erase(std::remove(uHits_cl.begin(),uHits_cl.end(),*h),uHits_cl.end());
	  }
	}
	if(FinalLine.size()==0 && LineHits.size()>0){
	    for(CP::THitSelection::iterator h=LineHits.begin();h!=LineHits.end();++h){
	    uHits_cl.erase(std::remove(uHits_cl.begin(),uHits_cl.end(),*h),uHits_cl.end());
	  }
	  }
  if(FinalLine.size()==0 && LineHits.size()<3){failHough=true;}
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

 if(uclusters->size()>0){
    for(CP::TReconObjectContainer::iterator it = uclusters->begin();it!=uclusters->end();++it){
      CP::THandle<CP::THitSelection> hits = (*it)->GetHits();
      if(hits){
	double avgDist=0;
	for(CP::THitSelection::iterator h = hits->begin();h!=(hits->end()-1);++h){avgDist+= ((*h)->GetPosition()- (*(h+1))->GetPosition()).Mag();
	}
	avgDist/=((double)hits->size()-1);
	distance2D->Fill(avgDist,(int)hits->size());
      }
    }
      }

      
      std::cout<<"NClusters="<<uclusters->size()<<std::endl;
     
      CaptNamedLog("TClusterUHits",
		   "With " << uclusters->size() << " clusters_uh");
    
    }


     
  if(vHits.size()>0)
    {
      CP::THitSelection vHits_cl;

      // vHits_cl=ClusteredHits(vHits,false,5,25);
      //	std::copy(vHits.begin(),vHits.end(),back_inserter(vHits_cl));
      for(CP::THitSelection::iterator h=vHits.begin();h !=vHits.end();++h){
	double hw=CP::GeomId::Captain::GetWireNumber((*h)->GetConstituent()->GetGeomId());
	 if(hw==308 || hw==310) continue;
	   vHits_cl.push_back(*h);
	 
      }
      bool failHough=false;
      typedef CP::THoughTrans< CP::THandle<THit> > HoughAlgorithm;

      while(!failHough){
	HoughAlgorithm Hough(2*180,10000);

	
	std::vector<std::pair<double,double>> points;

	for(CP::THitSelection::iterator h=vHits_cl.begin();h!=vHits_cl.end();++h){
	  double sign = 1.0;
	  if((*h)->GetPosition().X()<0)sign=-1.0;
	  double x=sign*sqrt((*h)->GetPosition().X()*(*h)->GetPosition().X()+(*h)->GetPosition().Y()*(*h)->GetPosition().Y());//(sin(30*3.14159265/180.0)*sin(30*3.14159265/180.0));
	  points.push_back(std::make_pair(x,(*h)->GetPosition().Z()));
	}
	
	Hough.HoughTransform(points.begin(),points.end());

	
	std::pair<double,double> lineParam=Hough.GetLineParam(3);

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
	  //  distance->Fill(dist);
	  if(dist<4*unit::mm) LineHits.push_back(*h);
	}

	CP::THitSelection FinalLine;
	
	if(LineHits.size()>2){
	  std::sort(LineHits.begin(),LineHits.end(),CompX);

	  	  FinalLine=CP::TClustering2D::ClusteredHits(LineHits,true,1,40*unit::mm);
	  //	std::copy(LineHits.begin(),LineHits.end(),back_inserter(FinalLine));
if(FinalLine.size()>finallinecut){
	  CP::THandle<CP::TReconCluster> cluster = CreateCluster("cluster_vh",FinalLine.begin(),FinalLine.end());
	  if(cluster){
#ifndef WIRE_CLUSTERING
	     CP::THandle<CP::THitSelection> hits = cluster->GetHits();
	     	double avgDist=0;
      if(hits){

	for(CP::THitSelection::iterator h = hits->begin();h!=(hits->end()-1);++h){avgDist+= ((*h)->GetPosition()- (*(h+1))->GetPosition()).Mag();
	}
      }
	if(cluster->GetHits()->size()>2){
	    vclusters->push_back(cluster);
	    final->push_back(cluster);
	    }
#endif
	  }
 }
	}
	
	if(FinalLine.size()>0){
	  for(CP::THitSelection::iterator h=FinalLine.begin();h!=FinalLine.end();++h){
	    vHits_cl.erase(std::remove(vHits_cl.begin(),vHits_cl.end(),*h),vHits_cl.end());
	  }
	}
        if(FinalLine.size()==0 && LineHits.size()>0){
	    for(CP::THitSelection::iterator h=LineHits.begin();h!=LineHits.end();++h){
	    vHits_cl.erase(std::remove(vHits_cl.begin(),vHits_cl.end(),*h),vHits_cl.end());
	  }
	  }
  if(FinalLine.size()==0 && LineHits.size()<3){failHough=true;}
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


       if(vclusters->size()>0){
    for(CP::TReconObjectContainer::iterator it = vclusters->begin();it!=vclusters->end();++it){
      CP::THandle<CP::THitSelection> hits = (*it)->GetHits();
      if(hits){
	double avgDist=0;
	for(CP::THitSelection::iterator h = hits->begin();h!=(hits->end()-1);++h){avgDist+= ((*h)->GetPosition()- (*(h+1))->GetPosition()).Mag();
	}
	avgDist/=((double)hits->size()-1);
	distance2D->Fill(avgDist,(int)hits->size());
      }
    }
      }
      std::cout<<"NClusters="<<vclusters->size()<<std::endl;
	
     
      CaptNamedLog("TClusterVHits",
		   "With " << vclusters->size() << " clusters_vh");
    
    }
    
   std::unique_ptr<TH2F> HitsX(new TH2F("HitsForX","HitsForX",340,0,340,9600,0,9600));
  // std::unique_ptr<TH2F> HitsX(new TH2F("HitsForX","HitsForX",1000,-500,500,2000,-1000,1000));
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
	  double sign = 1.0;
	  if((*h)->GetPosition().X()<0)sign=-1.0;
	   double x=sign*sqrt((*h)->GetConstituent()->GetPosition().X()*(*h)->GetConstituent()->GetPosition().X()+(*h)->GetConstituent()->GetPosition().Y()*(*h)->GetConstituent()->GetPosition().Y());
	   // i = HitsX->GetXaxis()->FindFixBin(x);
	   // std::cout<<(*h)->GetConstituent()->GetTime()*1.6/1000<<std::endl;
	   // j = HitsX->GetYaxis()->FindFixBin((*h)->GetConstituent()->GetTime()*1.6/1000);
	  HitsX->SetBinContent(i,j,coll*10);
	  
	}

      }
      ++coll;
    }

#ifdef DEBUG_HISTOGRAMS
    HitsX->Draw("COLZ");
    gPad->Print("plots/XHits_cl.C");
#endif
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
    
#ifdef DEBUG_HISTOGRAMS
    HitsU->Draw("COLZ");
    gPad->Print("plots/UHits_cl.C");
#endif
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
 
#ifdef DEBUG_HISTOGRAMS
    HitsV->Draw("COLZ");
    gPad->Print("plots/VHits_cl.C");
#endif
  }    
   
#ifdef DEBUG_HISTOGRAMS
  distance->Draw();
  gPad->Print("plots/distance.C");
  distance2D->Draw("COLZ");
  gPad->Print("plots/distance2D.C");
#endif
  result->AddResultsContainer(xclusters.release());
  result->AddResultsContainer(uclusters.release());
  result->AddResultsContainer(vclusters.release());
  result->AddHits(hitsAll.release());
  result->AddResultsContainer(final.release());
  
  return result;
}
