#include "TPathFollow.hxx"
#include "HitUtilities.hxx"
#include "CompareReconObjects.hxx"
#include "CreateTrack.hxx"
#include "CreateCluster.hxx"

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

std::pair<double,double> lineParam(TVector3 midle,TVector3 hit){
  double signM = 1.0;
  if(midle.X()<0)signM=-1.0;
  double xM=signM*sqrt(midle.X()*midle.X()+midle.Y()*midle.Y());
  double yM=midle.Z();
  double sign = 1.0;
  if(hit.X()<0)sign=-1.0;
  double x=sign*sqrt(hit.X()*hit.X()+hit.Y()*hit.Y());
  double y = hit.Z();

  double a = (yM-y)/(xM-x);
  double b = y-a*x;
  return std::make_pair(a,b);
  
}

double GetDistPath(double x0, double y0, std::pair<double,double> line)
{
  double a = line.first;
  double b = line.second;
  double d;
  d=fabs(a*x0-y0+b)/sqrt(a*a+1);

  return d;
}

  int CheckObjectPlaneP(const CP::THandle<CP::TReconBase>& obj){
    
 
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

    return -1;
  }

CP::TReconObjectContainer PathReCluster(CP::TReconObjectContainer& clusters, CP::THitSelection& unclusteredHits, int dir){
  CP::TReconObjectContainer newClusters;

  std::sort(clusters.begin(),clusters.end(),[](const CP::THandle<CP::TReconCluster> lcl,const CP::THandle<CP::TReconCluster> rcl){
      double lch = lcl->GetHits()->size();
      double rch = rcl->GetHits()->size();
      return lch>rch;
    });




  for(CP::TReconObjectContainer::iterator c=clusters.begin(); c != clusters.end();){
    CP::THandle<CP::TReconCluster> clust = (*c);
    TVector3 midle = clust->GetPosition().Vect();

    bool path = true;
    bool increment = true;
    bool foundH=false;
    bool foundCl=false;
    CP::THandle<CP::TReconCluster> tempC(new CP::TReconCluster);
    tempC=clust;
    
    while(path){
      
      midle = tempC->GetPosition().Vect();
      CP::THandle<CP::THitSelection> tempCHits = tempC->GetHits();
      std::sort(tempCHits->begin(),tempCHits->end(),[](const CP::THandle<CP::THit> lh,const CP::THandle<CP::THit> rh){
	    double lhX = lh->GetPosition().X();
	    double rhX = rh->GetPosition().X();
      return lhX<rhX;
    });
      
      TVector3 tempMax;
      if(dir==0){
      tempMax = tempCHits->back()->GetPosition();
      double ht=(tempCHits->back()->GetConstituent()->GetTime()+1.6*unit::ms)/(500*unit::ns);
	  double hw=CP::GeomId::Captain::GetWireNumber(tempCHits->back()->GetGeomId());
	  
      }
      else if(dir==1){
	tempMax = tempCHits->front()->GetPosition();
	double ht=(tempCHits->front()->GetConstituent()->GetTime()+1.6*unit::ms)/(500*unit::ns);
	  double hw=CP::GeomId::Captain::GetWireNumber(tempCHits->front()->GetGeomId());
	 
      }
      //define parameters for how far one looks from las hit 35/5 or 40/6.3 pure guess
      double cubel=40;
      double cubew=6.3;
    double distH=cubew;
    double distCl=cubew;
    double lengthH=9999;
    double lengthCl=9999;
    if(clusters.size()>1){
      for(CP::TReconObjectContainer::iterator c2=c+1; c2 != clusters.end();++c2){
	CP::THandle<CP::TReconCluster> clust2 = (*c2);
	CP::THandle<CP::THitSelection> c2Hits = clust2->GetHits();
	for(CP::THitSelection::iterator c2h=c2Hits->begin();c2h!=c2Hits->end();++c2h){
	  TVector3 hit = (*c2h)->GetPosition();
	  double sign =1.0;
	  if(hit.X()<0)sign=-1.0;
	  double x=sign*sqrt(hit.X()*hit.X()+hit.Y()*hit.Y());
	  double y = hit.Z();
	  double dist = GetDistPath(x,y,lineParam(midle,tempMax)) ;
       
	  if((hit-tempMax).Mag()<lengthCl){distCl=dist;lengthCl=(hit-tempMax).Mag();}
	}
      }
    }
    if(unclusteredHits.size()>0){
      for(CP::THitSelection::iterator uh=unclusteredHits.begin();uh!=unclusteredHits.end();++uh){
	TVector3 hit = (*uh)->GetPosition();
	double sign =1.0;
	if(hit.X()<0)sign=-1.0;
	double x=sign*sqrt(hit.X()*hit.X()+hit.Y()*hit.Y());
	double y = hit.Z();
	double dist = GetDistPath(x,y,lineParam(midle,tempMax)) ;
	if((hit-tempMax).Mag()<lengthH){distH=dist;lengthH=(hit-tempMax).Mag();}
      }
    }
    
    if(distH<distCl && distH<cubew && lengthH<cubel){
      foundH=true;
            for(CP::THitSelection::iterator uh=unclusteredHits.begin();uh!=unclusteredHits.end();++uh){
	TVector3 hit = (*uh)->GetPosition();
	double sign =1.0;
	if(hit.X()<0)sign=-1.0;
	double x=sign*sqrt(hit.X()*hit.X()+hit.Y()*hit.Y());
	double y = hit.Z();
	double dist = GetDistPath(x,y,lineParam(midle,tempMax)) ;
	if(fabs(dist-distH)<0.1 && fabs((hit-tempMax).Mag()-lengthH)<0.1){
	  
	  tempCHits->push_back(*uh);
	  CP::THandle<CP::TReconCluster> cluster = CreateCluster("tempCl",tempCHits->begin(),tempCHits->end());
	  tempC = cluster;
	  
	  unclusteredHits.erase(uh);
	  break;
	}
      }
    }else{foundH=false;}
    
    if((distCl<distH && distCl<cubew && lengthCl<cubel)||(!foundH && distCl<cubew && lengthCl<cubel)){
      foundCl=true;
      for(CP::TReconObjectContainer::iterator c2=c+1; c2 != clusters.end();++c2){
	CP::THandle<CP::TReconCluster> clust2 = (*c2);
	CP::THandle<CP::THitSelection> c2Hits = clust2->GetHits();
	bool stop=false;

	for(CP::THitSelection::iterator c2h=c2Hits->begin();c2h!=c2Hits->end();++c2h){
	  
	  TVector3 hit = (*c2h)->GetPosition();
	  double sign =1.0;
	  if(hit.X()<0)sign=-1.0;
	  double x=sign*sqrt(hit.X()*hit.X()+hit.Y()*hit.Y());
	  double y = hit.Z();
	  double dist = GetDistPath(x,y,lineParam(midle,tempMax)) ;

	  if(fabs(dist-distCl)<0.1 && fabs((hit-tempMax).Mag()-lengthCl)<0.1){
	    stop=true;
	    break;
	  }
	}
	if(stop){
	  CP::THitSelection remainingHits;
	  for(CP::THitSelection::iterator c2h=c2Hits->begin();c2h!=c2Hits->end();++c2h){
	    TVector3 hit = (*c2h)->GetPosition();
	  double sign =1.0;
	  if(hit.X()<0)sign=-1.0;
	  double x=sign*sqrt(hit.X()*hit.X()+hit.Y()*hit.Y());
	  double y = hit.Z();
	  double dist = GetDistPath(x,y,lineParam(midle,tempMax));
	  if(dist<cubew && (hit-tempMax).Mag()<cubel){
	    tempCHits->push_back(*c2h);
	  }else{remainingHits.push_back(*c2h);}
	  }

	  CP::THandle<CP::TReconCluster> clusterOld(new CP::TReconCluster);
	  bool leftCL=false;
	  if(remainingHits.size()>2){
	    clusterOld = CreateCluster("redoOld",remainingHits.begin(),remainingHits.end());
	    leftCL=true;
	  }else if(remainingHits.size()==1){
	     unclusteredHits.push_back(remainingHits[0]);
	  }else if(remainingHits.size()==2){
	    unclusteredHits.push_back(remainingHits[0]);
	    unclusteredHits.push_back(remainingHits[1]);
	  }
	 
	  CP::THandle<CP::TReconCluster> cluster = CreateCluster("tempCl",tempCHits->begin(),tempCHits->end());
	 
	  tempC = cluster;
	
	  if(leftCL){
	    (*c2)=clusterOld;
	  }else{
	   
	    clusters.erase(c2);
	  }

	  break;
	}
      }
    }else{foundCl=false;}
    
    if(!foundH && !foundCl)path=false;
    
    }

    if(tempC){newClusters.push_back(tempC);
      //  std::cout<<"newone"<<std::endl;
    }
    else{newClusters.push_back(*c);
      //std::cout<<"oldone"<<std::endl;
    }
    
    ++c;
  }

  
  return newClusters;

}



CP::TPathFollow::TPathFollow()
    : TAlgorithm("TPathFollow", 
                 "Path Follow algorithm") {

 
}

CP::TPathFollow::~TPathFollow() { }



CP::THandle<CP::TAlgorithmResult>
CP::TPathFollow::Process(const CP::TAlgorithmResult& input,
                               const CP::TAlgorithmResult& input1,
                               const CP::TAlgorithmResult&) {

  //ReconObjects 
  CP::THandle<CP::TReconObjectContainer> inputObjects 
      = input.GetResultsContainer();

  //allHits

  CP::THandle<CP::THitSelection> allDriftHits = input.GetHits();
  
  //1
    CaptLog("TPathFollow Process " << GetEvent().GetContext());

    if (!inputObjects) {
        CaptError("No input objects");
        return CP::THandle<CP::TAlgorithmResult>();
    }
    
    // Create the output containers.
    CP::THandle<CP::TAlgorithmResult> result = CreateResult();
    std::unique_ptr<CP::TReconObjectContainer> 
        final(new CP::TReconObjectContainer("final"));
        std::unique_ptr<CP::TReconObjectContainer> 
        clustersX_pf(new CP::TReconObjectContainer("clustersX_pf"));
	    std::unique_ptr<CP::TReconObjectContainer> 
        clustersU_pf(new CP::TReconObjectContainer("clustersU_pf"));
	    std::unique_ptr<CP::TReconObjectContainer> 
        clustersV_pf(new CP::TReconObjectContainer("clustersV_pf"));

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
	int plane = CheckObjectPlaneP(cluster);
	if(plane==1){
	  clustersX.push_back(*cl);}
	else if(plane==2){
	  clustersU.push_back(*cl);}
	else if(plane==3){
	  clustersV.push_back(*cl);}
	else{
            std::cout<<"PLANEDEFININGFORCLUSTERSDOESNOTWORK"<<std::endl;
        }
    }
    
	// CP::THandle<CP::THitSelection> allDriftHits = GetEvent().Get<CP::THitSelection>("~/fits/TCaptainRecon/THitTransfer/final");
    if(clustersX.size()>0)
      {
	CP::THitSelection unclusteredHits;
	std::set<CP::THandle<CP::THit>> allclustXHits;
	for(CP::TReconObjectContainer::iterator cl = clustersX.begin(); cl!= clustersX.end();++cl)
	  { CP::THandle<CP::TReconCluster> cluster = *cl;
	    CP::THandle<CP::THitSelection> clHits = cluster->GetHits();
	    for(CP::THitSelection::iterator h = clHits->begin();h!=clHits->end();++h){
	      allclustXHits.insert(*h);
	    }
	  }
	      for(CP::THitSelection::iterator ha = allDriftHits->begin();ha!=allDriftHits->end();++ha){
		if(CP::GeomId::Captain::GetWirePlane((*ha)->GetGeomId())==CP::GeomId::Captain::kXPlane){
		   if (allclustXHits.find(*ha) != allclustXHits.end()) continue;
		  unclusteredHits.push_back(*ha);
		}
	      }
	      
	      
              // std::cout<<"Start: clusters="<<clustersX.size()<<"; unusedHits="<<unclusteredHits.size()<<std::endl;
	clustersX=PathReCluster(clustersX,unclusteredHits,0);
	// std::cout<<"AfterForwardPath: clusters="<<clustersX.size()<<"; unusedHits="<<unclusteredHits.size()<<std::endl;
        clustersX=PathReCluster(clustersX,unclusteredHits,1);
	// std::cout<<"AfterBackwardPath: clusters="<<clustersX.size()<<"; unusedHits="<<unclusteredHits.size()<<std::endl;
	for(CP::TReconObjectContainer::iterator cl = clustersX.begin();cl!=clustersX.end();++cl){
	   CP::THandle<CP::THitSelection> clHits = (*cl)->GetHits();
	 std::sort(clHits->begin(),clHits->end(),[](const CP::THandle<CP::THit> lh,const CP::THandle<CP::THit> rh){
	    double lhX = lh->GetPosition().X();
	    double rhX = rh->GetPosition().X();
      return lhX<rhX;
    });
	 CP::THitSelection::iterator h ;
	 h=clHits->begin();
	 double cutDist = ((*h)->GetPosition()-(*(h+1))->GetPosition()).Mag();
	 if(cutDist>35)clHits->erase(h);
	 h=clHits->end();
	 cutDist = ((*(h-1))->GetPosition()-(*(h-2))->GetPosition()).Mag();
	 if(cutDist>35)clHits->erase(h-1);

	 CP::THandle<CP::TReconCluster> clFinal = CreateCluster("finalclusterX_pf",clHits->begin(),clHits->end());

	 if(clFinal){
	   if(clFinal->GetHits()->size()>3){
	  clustersX_pf->push_back(clFinal);
	  final->push_back(clFinal);
	   }
	 }
	}
      }

    CaptNamedLog("TPathFollowXClusters",
                 "With " << (*clustersX_pf).size() << " clustersX_pf"
                 << " from " << clustersX.size() << " x clusters");

     if(clustersU.size()>0)
      {
        
	CP::THitSelection unclusteredHits;
	std::set<CP::THandle<CP::THit>> allclustUHits;
	for(CP::TReconObjectContainer::iterator cl = clustersU.begin(); cl!= clustersU.end();++cl)
	  { CP::THandle<CP::TReconCluster> cluster = *cl;
	    CP::THandle<CP::THitSelection> clHits = cluster->GetHits();
	    for(CP::THitSelection::iterator h = clHits->begin();h!=clHits->end();++h){
	      allclustUHits.insert(*h);
	    }
	  }
	      for(CP::THitSelection::iterator ha = allDriftHits->begin();ha!=allDriftHits->end();++ha){
		if(CP::GeomId::Captain::GetWirePlane((*ha)->GetGeomId())==CP::GeomId::Captain::kUPlane){
		   if (allclustUHits.find(*ha) != allclustUHits.end()) continue;
		  unclusteredHits.push_back(*ha);
		}
	      }
              // std::cout<<"Start: clustersU="<<clustersU.size()<<"; unusedHits="<<unclusteredHits.size()<<std::endl;
	clustersU=PathReCluster(clustersU,unclusteredHits,0);
	// std::cout<<"AfterForwardPath: clustersU="<<clustersU.size()<<"; unusedHits="<<unclusteredHits.size()<<std::endl;
        clustersU=PathReCluster(clustersU,unclusteredHits,1);
	// std::cout<<"AfterBackwardPath: clusters="<<clustersU.size()<<"; unusedHits="<<unclusteredHits.size()<<std::endl;
	for(CP::TReconObjectContainer::iterator cl = clustersU.begin();cl!=clustersU.end();++cl){
	  
	   CP::THandle<CP::THitSelection> clHits = (*cl)->GetHits();
	  
	 std::sort(clHits->begin(),clHits->end(),[](const CP::THandle<CP::THit> lh,const CP::THandle<CP::THit> rh){
	    double lhX = lh->GetPosition().X();
	    double rhX = rh->GetPosition().X();
      return lhX<rhX;
    });
	 CP::THitSelection::iterator h ;
	 h=clHits->begin();
	 double cutDist = ((*h)->GetPosition()-(*(h+1))->GetPosition()).Mag();
	 if(cutDist>35)clHits->erase(h);
	 h=clHits->end();
	 cutDist = ((*(h-1))->GetPosition()-(*(h-2))->GetPosition()).Mag();
	 if(cutDist>35)clHits->erase(h-1);

	 CP::THandle<CP::TReconCluster> clFinal = CreateCluster("finalclusterU_pf",clHits->begin(),clHits->end());

	 if(clFinal){
	   if(clFinal->GetHits()->size()>3){
	  clustersU_pf->push_back(clFinal);
	  final->push_back(clFinal);
	   }
	 }
	
	}
      }

      CaptNamedLog("TPathFollowUClusters",
                 "With " << (*clustersU_pf).size() << " clustersU_pf"
                 << " from " << clustersU.size() << " u clusters");
      
    if(clustersV.size()>0)
      {
        
	CP::THitSelection unclusteredHits;
	std::set<CP::THandle<CP::THit>> allclustVHits;
	for(CP::TReconObjectContainer::iterator cl = clustersV.begin(); cl!= clustersV.end();++cl)
	  { CP::THandle<CP::TReconCluster> cluster = *cl;
	    CP::THandle<CP::THitSelection> clHits = cluster->GetHits();
	    for(CP::THitSelection::iterator h = clHits->begin();h!=clHits->end();++h){
	      allclustVHits.insert(*h);
	    }
	  }
	      for(CP::THitSelection::iterator ha = allDriftHits->begin();ha!=allDriftHits->end();++ha){
		if(CP::GeomId::Captain::GetWirePlane((*ha)->GetGeomId())==CP::GeomId::Captain::kVPlane){
		   if (allclustVHits.find(*ha) != allclustVHits.end()) continue;
		  unclusteredHits.push_back(*ha);
		}
	      }
	  
              // std::cout<<"Start: clustersV="<<clustersV.size()<<"; unusedHits="<<unclusteredHits.size()<<std::endl;
	clustersV=PathReCluster(clustersV,unclusteredHits,0);
	// std::cout<<"AfterForwardPath: clustersV="<<clustersV.size()<<"; unusedHits="<<unclusteredHits.size()<<std::endl;
        clustersV=PathReCluster(clustersV,unclusteredHits,1);
	// std::cout<<"AfterBackwardPath: clustersV="<<clustersV.size()<<"; unusedHits="<<unclusteredHits.size()<<std::endl;
	for(CP::TReconObjectContainer::iterator cl = clustersV.begin();cl!=clustersV.end();++cl){
	  CP::THandle<CP::THitSelection> clHits = (*cl)->GetHits();
	  // std::cout<<clHits->size()<<std::endl;
	 std::sort(clHits->begin(),clHits->end(),[](const CP::THandle<CP::THit> lh,const CP::THandle<CP::THit> rh){
	    double lhX = lh->GetPosition().X();
	    double rhX = rh->GetPosition().X();
      return lhX<rhX;
    });
	 CP::THitSelection::iterator h ;
	 h=clHits->begin();
	 double cutDist = ((*h)->GetPosition()-(*(h+1))->GetPosition()).Mag();
	 if(cutDist>35)clHits->erase(h);
	 h=clHits->end();
	 cutDist = ((*(h-1))->GetPosition()-(*(h-2))->GetPosition()).Mag();
	 if(cutDist>35)clHits->erase(h-1);
	 
	 CP::THandle<CP::TReconCluster> clFinal = CreateCluster("finalclusterV_pf",clHits->begin(),clHits->end());

	 if(clFinal){
	   if(clFinal->GetHits()->size()>3){
	  clustersV_pf->push_back(clFinal);
	  final->push_back(clFinal);
	   }
	 }
	}
      }

     CaptNamedLog("TPathFollowVClusters",
                 "With " << (*clustersV_pf).size() << " clustersV_pf"
                 << " from " << clustersV.size() << " v clusters");

     std::unique_ptr<TH2F> HitsX(new TH2F("HitsForX","HitsForX",340,0,340,9600,0,9600));
  std::unique_ptr<TH2F> HitsU(new TH2F("HitsForU","HitsForU",340,0,340,9600,0,9600));
  std::unique_ptr<TH2F>  HitsV(new TH2F("HitsForV","HitsForV",340,0,340,9600,0,9600));
  
  if(clustersX_pf->size()>0){
    int coll=1;
    for(CP::TReconObjectContainer::iterator it = clustersX_pf->begin();it!=clustersX_pf->end();++it){
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
    gPad->Print("plots/XHits_pf.C");
#endif
  }

  if(clustersU_pf->size()>0){
    int coll=1;
    for(CP::TReconObjectContainer::iterator it = clustersU_pf->begin();it!=clustersU_pf->end();++it){
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
	  HitsU->SetBinContent(i,j,coll*10);
	  
	}

      }
      ++coll;
    }

#ifdef DEBUG_HISTOGRAMS
    HitsU->Draw("COLZ");
    gPad->Print("plots/UHits_pf.C");
#endif
  }

  if(clustersV_pf->size()>0){
    int coll=1;
    for(CP::TReconObjectContainer::iterator it = clustersV_pf->begin();it!=clustersV_pf->end();++it){
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
	  HitsV->SetBinContent(i,j,coll*10);
	  
	}

      }
      ++coll;
    }

#ifdef DEBUG_HISTOGRAMS
    HitsV->Draw("COLZ");
    gPad->Print("plots/VHits_pf.C");
#endif
  }


    result->AddResultsContainer(final.release());

    return result;
}
