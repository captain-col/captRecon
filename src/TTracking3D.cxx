#include "TTracking3D.hxx"
#include "HitUtilities.hxx"
#include "CompareReconObjects.hxx"
#include "CreateTrack.hxx"
#include "CreateHit.hxx"
#include "HitConnector.hxx"

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
#include <string.h>
#include <sstream>

#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"

 std::string toString(int i)
  {
    std::ostringstream s;
    s << i;
    return s.str();
  }

  int CheckObjectPlaneT(const CP::THandle<CP::TReconTrack>& obj){
    
 
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


bool CompHitsW(const CP::THandle<CP::THit>& lhs, const CP::THandle<CP::THit>& rhs){
return CP::GeomId::Captain::GetWireNumber(lhs->GetGeomId()) < CP::GeomId::Captain::GetWireNumber(rhs->GetGeomId());
}

bool CompHitsWB(const CP::THandle<CP::THit>& lhs, const CP::THandle<CP::THit>& rhs){
return CP::GeomId::Captain::GetWireNumber(lhs->GetGeomId()) > CP::GeomId::Captain::GetWireNumber(rhs->GetGeomId());
}

bool CompHitsZ(const CP::THandle<CP::THit>& lhs, const CP::THandle<CP::THit>& rhs){
  return lhs->GetPosition().Z() < rhs->GetPosition().Z();
}

bool CompHitsT(const CP::THandle<CP::THit>& lhs, const CP::THandle<CP::THit>& rhs){
  return lhs->GetTime() < rhs->GetTime();
}

bool CheckPosition(TVector3 s, TVector3 f, TVector3 posX){
  if(posX.Z()>s.Z() && posX.Z()<f.Z()) return true;
   if(posX.Z()<s.Z() && posX.Z()>f.Z()) return true;
  return false;
}
int WireDir( CP::THandle<CP::TReconTrack> trackX, CP::THandle<CP::TReconTrack> trackU,CP::THandle<CP::TReconTrack> trackV){
  
  CP::THandle<CP::THitSelection> hitX = trackX->GetHits();
  CP::THandle<CP::THitSelection> hitU = trackU->GetHits();
  std::set<int> nWire;
    for(CP::THitSelection::iterator h = (*hitX).begin() ; h!=(*hitX).end();++h)
    {
      nWire.insert(CP::GeomId::Captain::GetWireNumber((*h)->GetGeomId()));
    }
int nWireX=nWire.size();
  nWire.clear();

   for(CP::THitSelection::iterator h = (*hitU).begin() ; h!=(*hitU).end();++h)
    {
      nWire.insert(CP::GeomId::Captain::GetWireNumber((*h)->GetGeomId()));
    }
 int nWireU=nWire.size();
  nWire.clear();

  if(nWireX<3 || nWireU<3)return 10;

    std::sort((*hitX).begin(),(*hitX).end(),CompHitsW);
 
  std::sort((*hitU).begin(),(*hitU).end(),CompHitsW);
  
  TVector3 xuFF = PositionXY((*hitX).front()->GetConstituent(),(*hitU).front()->GetConstituent());
      TVector3 xuFB = PositionXY((*hitX).front()->GetConstituent(),(*hitU).back()->GetConstituent());

      if(trackV){    
  CP::THandle<CP::THitSelection> hitV = trackV->GetHits();

   for(CP::THitSelection::iterator h = (*hitV).begin() ; h!=(*hitV).end();++h)
    {
      nWire.insert(CP::GeomId::Captain::GetWireNumber((*h)->GetGeomId()));
    }
 int nWireV=nWire.size();
  nWire.clear();

  if(nWireV<3)return 10;
  
  std::sort((*hitV).begin(),(*hitV).end(),CompHitsW);
  TVector3 xvFF = PositionXY((*hitX).front()->GetConstituent(),(*hitV).front()->GetConstituent());
  TVector3 xvFB = PositionXY((*hitX).front()->GetConstituent(),(*hitV).back()->GetConstituent());
 
  double dist1 = (xuFF-xvFF).Mag();
  double dist2 = (xuFF-xvFB).Mag();
  double dist3 = (xuFB-xvFB).Mag();
  double dist4 = (xuFB-xvFF).Mag();
  double min = std::min(dist1,dist2);
  min = std::min(min,dist3);
  min = std::min(min,dist4);

  if(min==dist1){return 0;}
  if(min==dist2){return 1;}
  if(min==dist3){return 2;}
  if(min==dist4){return 3;}
      }else{
	double dist1 = (xuFF-(*hitX).front()->GetPosition()).Mag();
	double dist2 = (xuFB-(*hitX).front()->GetPosition()).Mag();
	if(dist1>dist2 || dist1==dist2){return 0;}
	if(dist1<dist2){return 1;}
      }
  
  return 9999;
  }

double MaxZ(CP::THandle<CP::TReconTrack> track){
  CP::THandle<CP::THitSelection> hits = track->GetHits();
  double maxZ=-99999;
  for (CP::THitSelection::iterator h = (*hits).begin();
       h != (*hits).end(); ++h) {
    double z = (*h)->GetPosition().Z();
    maxZ=std::max(maxZ,z);
  }
  return maxZ;
}

double MinZ(CP::THandle<CP::TReconTrack> track){
  CP::THandle<CP::THitSelection> hits = track->GetHits();
  double minZ=99999;
  for (CP::THitSelection::iterator h = (*hits).begin();
       h != (*hits).end(); ++h) {
    double z = (*h)->GetPosition().Z();
    minZ=std::min(minZ,z);
  }
  return minZ;
}

struct CompTracks{

  CompTracks(double maxZX, double minZX){this->maxZX=maxZX;this->minZX=minZX;}

  bool operator()(const CP::THandle<CP::TReconTrack>& lhs, const CP::THandle<CP::TReconTrack>& rhs)
	      {
		double l=abs(MaxZ(lhs)-maxZX)+abs(MinZ(lhs)-minZX);
		double r=abs(MaxZ(rhs)-maxZX)+abs(MinZ(rhs)-minZX);
		return l<r;}
  
  double maxZX;
  double minZX;
};

bool Assemble3DTrack( CP::THandle<CP::TReconTrack> trackX, CP::THandle<CP::TReconTrack> trackU, CP::THandle<CP::TReconTrack> trackV,CP::TReconObjectContainer& match3,int trackNum){

  
  CP::THandle<CP::THitSelection> hitX_t = trackX->GetHits();
  CP::THandle<CP::THitSelection> hitU_t = trackU->GetHits();
  CP::THandle<CP::THitSelection> hitV_t = trackV->GetHits();
  std::set<CP::THandle<CP::THit>> setX;
  std::set<CP::THandle<CP::THit>> setU;
  std::set<CP::THandle<CP::THit>> setV;
  CP::THandle<CP::THitSelection> hitX(new CP::THitSelection);
  CP::THandle<CP::THitSelection> hitU(new CP::THitSelection);
  CP::THandle<CP::THitSelection> hitV(new CP::THitSelection);
  TH2F* HitsX = new TH2F("HitsForX","HitsForX",340,0,340,9600,0,9600);
  TH2F* HitsU = new TH2F("HitsForU","HitsForU",340,0,340,9600,0,9600);
  TH2F* HitsV = new TH2F("HitsForV","HitsForV",340,0,340,9600,0,9600);

  
  //********************************************************************************
  //Since 3D hit finding algoritm later on operates with 2D hits, we extract unique constituents from "pseudo" 3D
  // Hits, formed for each 2D track
  //********************************************************************************

  
  for(CP::THitSelection::iterator h = (*hitX_t).begin() ; h!=(*hitX_t).end();++h){
  
   double ht=((*h)->GetConstituent()->GetTime()+1.6*unit::ms)/(500*unit::ns);
   double hw=CP::GeomId::Captain::GetWireNumber((*h)->GetGeomId());
   HitsX->Fill(hw,ht);
   CP::THandle<CP::THit> ch=(*h)->GetConstituent();
    setX.insert(ch);
    }
 HitsX->Draw();
 std::string nameX="plots/XHits_3D_"+toString(trackNum)+".C";
   gPad->Print(nameX.c_str());
  for(CP::THitSelection::iterator h = (*hitU_t).begin() ; h!=(*hitU_t).end();++h){
        double ht=((*h)->GetConstituent()->GetTime()+1.6*unit::ms)/(500*unit::ns);
	    double hw=CP::GeomId::Captain::GetWireNumber((*h)->GetGeomId());
	    HitsU->Fill(hw,ht);
     CP::THandle<CP::THit> ch=(*h)->GetConstituent();
      setU.insert(ch);
    }
 HitsU->Draw();
 std::string nameU="plots/UHits_3D_"+toString(trackNum)+".C";
   gPad->Print(nameU.c_str());
  for(CP::THitSelection::iterator h = (*hitV_t).begin() ; h!=(*hitV_t).end();++h){
        double ht=((*h)->GetConstituent()->GetTime()+1.6*unit::ms)/(500*unit::ns);
	    double hw=CP::GeomId::Captain::GetWireNumber((*h)->GetGeomId());
	    HitsV->Fill(hw,ht);
      CP::THandle<CP::THit> ch=(*h)->GetConstituent();
      setV.insert(ch);
    }
   HitsV->Draw();
   std::string nameV="plots/VHits_3D_"+toString(trackNum)+".C";
   gPad->Print(nameV.c_str());

  for(std::set<CP::THandle<CP::THit>>::iterator it = setX.begin();it!=setX.end();++it){
 
      hitX->push_back(*it); 
    }
   for(std::set<CP::THandle<CP::THit>>::iterator it = setU.begin();it!=setU.end();++it){
   
      hitU->push_back(*it);
    }
    for(std::set<CP::THandle<CP::THit>>::iterator it = setV.begin();it!=setV.end();++it){
    
      hitV->push_back(*it);
    }
  


  
  std::cout<<"#XHits="<<hitX->size()<<"; #UHis="<<hitU->size()<<"; #VHits="<<hitV->size()<<std::endl;


    //********************************************************************************
  //Hits should be sorted in a right order, so the geometric direction is determinde here(look at WireDir function   //for more deteils)
  //********************************************************************************
  
  int dir=WireDir(trackX,trackU,trackV);


  std::cout<<"dir="<<dir<<std::endl;
  std::sort((*hitX).begin(),(*hitX).end(),CompHitsW);
  if(dir==0){
    std::sort((*hitU).begin(),(*hitU).end(),CompHitsW);
    std::sort((*hitV).begin(),(*hitV).end(),CompHitsW);
  }
  if(dir==1){
    std::sort((*hitU).begin(),(*hitU).end(),CompHitsW);
    std::sort((*hitV).begin(),(*hitV).end(),CompHitsWB);
  }
  if(dir==2){
    std::sort((*hitU).begin(),(*hitU).end(),CompHitsWB);
    std::sort((*hitV).begin(),(*hitV).end(),CompHitsWB);
  }
  if(dir==3){
    std::sort((*hitU).begin(),(*hitU).end(),CompHitsWB);
    std::sort((*hitV).begin(),(*hitV).end(),CompHitsW);
  }
  if(dir==10){
    std::sort((*hitX).begin(),(*hitX).end(),CompHitsT);
     std::sort((*hitU).begin(),(*hitU).end(),CompHitsT);
    std::sort((*hitV).begin(),(*hitV).end(),CompHitsT);
  }
  
  //********************************************************************************
  //In order to determin how many parts should be in a track we find number of unique wires in each track
  // and set it as a baseline 
  //********************************************************************************
  int nWireX=0;
  int nWireU=0;
  int nWireV=0;
  std::set<int> nWire;
  for(CP::THitSelection::iterator h = (*hitX).begin() ; h!=(*hitX).end();++h)
    {
      nWire.insert(CP::GeomId::Captain::GetWireNumber((*h)->GetGeomId()));
    }
  nWireX=nWire.size();
  nWire.clear();
  for(CP::THitSelection::iterator h = (*hitU).begin() ; h!=(*hitU).end();++h)
    {
      nWire.insert(CP::GeomId::Captain::GetWireNumber((*h)->GetGeomId()));
    }
  nWireU=nWire.size();
  nWire.clear();
  for(CP::THitSelection::iterator h = (*hitV).begin() ; h!=(*hitV).end();++h)
    {
      nWire.insert(CP::GeomId::Captain::GetWireNumber((*h)->GetGeomId()));
    }
  nWireV=nWire.size();
  nWire.clear();

  int minSepar = std::min(nWireX,nWireU);
   minSepar = std::min(nWireV,minSepar);


  if(minSepar>3) minSepar/=2;
  if(minSepar > 20 && minSepar <= 50){minSepar/=2;}
  if(minSepar > 50 && minSepar <= 100){minSepar/=5;}
  if(minSepar > 100){minSepar/=10;}
  int nparts=minSepar;
      CP::TReconObjectContainer writableClusters;
 


 int counterX=(int)(*hitX).size()/minSepar;
 int nMissedX = (*hitX).size()-counterX*minSepar;
 int counterU=(int)(*hitU).size()/minSepar;
  int nMissedU = (*hitU).size()-counterU*minSepar;
 int counterV=(int)(*hitV).size()/minSepar;
  int nMissedV = (*hitV).size()-counterV*minSepar;

    //********************************************************************************
  //Main algorithm: it takes parts of tracks and apply HitConnector to them / optionally it can just create all possible hits from given sets
  //********************************************************************************
	for(std::size_t l=0;l<nparts;++l){

    CP::THitSelection writableHits;
 CP::THitSelection inHitsX;
 CP::THitSelection inHitsU;
 CP::THitSelection inHitsV;

 CP::THandle<CP::THitSelection> inHitsXH(new CP::THitSelection);
 CP::THandle<CP::THitSelection> inHitsUH(new CP::THitSelection);
 CP::THandle<CP::THitSelection> inHitsVH(new CP::THitSelection);
 // std::cout<<"l part="<<l<<std::endl;

    for(std::size_t j=l*counterX; j<counterX*(l+1) ; ++j){
      //std::cout<<"j part="<<j<<std::endl;
      inHitsX.push_back((*hitX)[j]);
      inHitsXH->push_back((*hitX)[j]);
    }
    for(std::size_t j=l*counterU; j<counterU*(l+1) ; ++j){
      inHitsU.push_back((*hitU)[j]);
      inHitsUH->push_back((*hitU)[j]);
    }
    for(std::size_t j=l*counterV; j<counterV*(l+1) ; ++j){
      inHitsV.push_back((*hitV)[j]);
      inHitsVH->push_back((*hitV)[j]);
    }
    if(l==(nparts-1)){
      if(nMissedX>0){
	int bound = (*hitX).size()-nMissedX;
	for(int j=bound; j<(*hitX).size() ; ++j){
	  inHitsX.push_back((*hitX)[j]);
	  inHitsXH->push_back((*hitX)[j]);
	}
      }
            if(nMissedU>0){
	int bound = (*hitU).size()-nMissedU;
	for(int j=bound; j<(*hitU).size() ; ++j){
	  inHitsU.push_back((*hitU)[j]);
	  inHitsUH->push_back((*hitU)[j]);
	}
      }
	          if(nMissedV>0){
	int bound = (*hitV).size()-nMissedV;
	for(int j=bound; j<(*hitV).size() ; ++j){
	  inHitsV.push_back((*hitV)[j]);
	  inHitsVH->push_back((*hitV)[j]);
	}
      }     
    }

 #define Hit_Con_Alg
 #ifdef Hit_Con_Alg
  writableHits = HitConnector3D(inHitsXH,inHitsUH,inHitsVH);
  #endif
  
 #ifdef All_3D_Hits
  if(inHitsX.size()>0){
   for(std::size_t i=0;i<inHitsX.size();++i){
 
     if(inHitsU.size()>0 ){
     for(std::size_t j=0;j<inHitsU.size();++j){

       if(inHitsV.size()>0 ){
	 for(std::size_t k=0;k<inHitsV.size();++k){

	   CP::THandle<CP::THit> h1=inHitsX[i];
	   CP::THandle<CP::THit> h2=inHitsU[j];
	   CP::THandle<CP::THit> h3=inHitsV[k];
	   TVector3 hitPosition;
	   OverlapedXY(h1,h2,h3,hitPosition);
	   if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
	    
	 }
       }else{
       CP::THandle<CP::THit> h1=inHitsX[i];
       CP::THandle<CP::THit> h2=inHitsU[j];
       if(!CreateHit(writableHits,PositionXY(h1,h2),0.0,h1,h2,CP::THandle<CP::THit>()))continue;}
     }
   }
   if(inHitsV.size()>0 && inHitsU.size()==0){
     for(std::size_t k=0;k<inHitsV.size();++k){
   
       CP::THandle<CP::THit> h1=inHitsX[i];
       CP::THandle<CP::THit> h2=inHitsV[k];
       if(!CreateHit(writableHits,PositionXY(h1,h2),0.0,h1,h2,CP::THandle<CP::THit>()))continue;
     }
   }
   }
 }
  #endif
  
 std::cout<<"#3DHits="<<writableHits.size()<<std::endl;
 
 if(writableHits.size()>0){


  //********************************************************************************
  //THis part might not be nessesary, but some hits might still have nan for time and Z, so put it here for now
  //********************************************************************************
 for(CP::THitSelection::iterator wh=writableHits.begin();wh!=writableHits.end();)
     {
     
       
       if(std::isnan((*wh)->GetPosition().Z()) || std::isnan(-1*(*wh)->GetPosition().Z()) || (*wh)->GetUncertainty().Z()==0){writableHits.erase(wh);}
       else++wh;
       
       }
   
   if(writableHits.size()>0){
   CP::THandle<CP::TReconCluster> wcl = CreateCluster("fortrack",writableHits.begin(),writableHits.end());
   writableClusters.push_back(wcl);

   }
 }
  writableHits.clear();
  }

  //********************************************************************************
  //From CLusters, formed in previos algorithm we form  track
  //********************************************************************************

	 if(writableClusters.size()>0){
   
    CP::THandle<CP::TReconTrack> track(new CP::TReconTrack);
     track= CreateTrackFromClusters("TTracking3D",writableClusters.begin(), writableClusters.end());

    if(track){
    match3.push_back(track);
    }

    return true;
  }

  return false;
  
}

bool Assemble2DTrack( CP::THandle<CP::TReconTrack> trackX, CP::THandle<CP::TReconTrack> trackU,CP::TReconObjectContainer& match2,int trackNum){

   //********************************************************************************
  //Algorithm follows the same logic as Assemble 3D, so check the other one for more detailed comments
  //********************************************************************************
 
  CP::THandle<CP::THitSelection> hitX_t = trackX->GetHits();
  CP::THandle<CP::THitSelection> hitU_t = trackU->GetHits();
  std::set<CP::THandle<CP::THit>> setX;
  std::set<CP::THandle<CP::THit>> setU;
  CP::THandle<CP::THitSelection> hitX(new CP::THitSelection);
  CP::THandle<CP::THitSelection> hitU(new CP::THitSelection);

   TH2F* HitsX = new TH2F("HitsForX","HitsForX",340,0,340,9600,0,9600);
    TH2F* HitsU = new TH2F("HitsForU","HitsForU",340,0,340,9600,0,9600);
  for(CP::THitSelection::iterator h = (*hitX_t).begin() ; h!=(*hitX_t).end();++h)
    {

       double ht=((*h)->GetConstituent()->GetTime()+1.6*unit::ms)/(500*unit::ns);
	    double hw=CP::GeomId::Captain::GetWireNumber((*h)->GetGeomId());
	    HitsX->Fill(hw,ht);
      CP::THandle<CP::THit> ch=(*h)->GetConstituent();
      setX.insert(ch);
    }

   HitsX->Draw();
   std::string nameX="plots/XHits_2D_"+toString(trackNum)+".C";
   gPad->Print(nameX.c_str());

  for(CP::THitSelection::iterator h = (*hitU_t).begin() ; h!=(*hitU_t).end();++h)
    {
       double ht=((*h)->GetConstituent()->GetTime()+1.6*unit::ms)/(500*unit::ns);
	    double hw=CP::GeomId::Captain::GetWireNumber((*h)->GetGeomId());
	    HitsU->Fill(hw,ht);
      CP::THandle<CP::THit> ch=(*h)->GetConstituent();
      setU.insert(ch);
    }

   HitsU->Draw();
   std::string nameU="plots/UHits_2D_"+toString(trackNum)+".C";
   gPad->Print(nameU.c_str());

 
  for(std::set<CP::THandle<CP::THit>>::iterator it = setX.begin();it!=setX.end();++it)
    {
      hitX->push_back(*it); 
    }
   for(std::set<CP::THandle<CP::THit>>::iterator it = setU.begin();it!=setU.end();++it)
    {
      hitU->push_back(*it);
    }


  std::sort((*hitX).begin(),(*hitX).end(),CompHitsT);
 
    std::sort((*hitU).begin(),(*hitU).end(),CompHitsT);

 
 int minSepar = std::min((int)(*hitX).size(),(int)(*hitU).size());

  if(minSepar>3) minSepar/=2;
  if(minSepar > 20 && minSepar <= 50){minSepar/=2;}
  if(minSepar > 50 && minSepar <= 100){minSepar/=5;}
  if(minSepar > 100){minSepar/=10;}
  int nparts=minSepar;
  
      CP::TReconObjectContainer writableClusters;



 int counterX=(int)(*hitX).size()/minSepar;
 int nMissedX = (*hitX).size()-counterX*minSepar;
 int counterU=(int)(*hitU).size()/minSepar;
  int nMissedU = (*hitU).size()-counterU*minSepar;


	for(std::size_t l=0;l<nparts;++l){

    CP::THitSelection writableHits;
 CP::THitSelection inHitsX;
 CP::THitSelection inHitsU;

 CP::THandle<CP::THitSelection> inHitsXH(new CP::THitSelection);
 CP::THandle<CP::THitSelection> inHitsUH(new CP::THitSelection);
 
 

    for(std::size_t j=l*counterX; j<counterX*(l+1) ; ++j){
   
      inHitsX.push_back((*hitX)[j]);
      inHitsXH->push_back((*hitX)[j]);
    }
    for(std::size_t j=l*counterU; j<counterU*(l+1) ; ++j){
      inHitsU.push_back((*hitU)[j]);
      inHitsUH->push_back((*hitU)[j]);
    }

    if(l==(nparts-1)){
      if(nMissedX>0){
	int bound = (*hitX).size()-nMissedX;
	for(int j=bound; j<(*hitX).size() ; ++j){
	  inHitsX.push_back((*hitX)[j]);
	  inHitsXH->push_back((*hitX)[j]);
	}
      }
            if(nMissedU>0){
	int bound = (*hitU).size()-nMissedU;
	for(int j=bound; j<(*hitU).size() ; ++j){
	  inHitsU.push_back((*hitU)[j]);
	  inHitsUH->push_back((*hitU)[j]);
	}
      }           
    }
    
    #define Hit_Con_Alg
#ifdef Hit_Con_Alg
    
 writableHits = HitConnector2D(inHitsXH,inHitsUH);

#endif


 //#define All_3D_Hits
#ifdef All_3D_Hits 
  if(inHitsX.size()>0){
   TVector3 hitPosition;
   for(std::size_t i=0;i<inHitsX.size();++i){
     
     if(inHitsU.size()>0){
     for(std::size_t k=0;k<inHitsU.size();++k){
       CP::THandle<CP::THit> h1=inHitsX[i];
       CP::THandle<CP::THit> h2=inHitsU[k];
       if(!CreateHit(writableHits,PositionXY(h1,h2),0.0,h1,h2,CP::THandle<CP::THit>()))continue;
     }
   }
   }
 }
#endif
 
 if(writableHits.size()>0){

    for(CP::THitSelection::iterator wh=writableHits.begin();wh!=writableHits.end();)
     {
       
       if(std::isnan((*wh)->GetPosition().Z()) || std::isnan(-1*(*wh)->GetPosition().Z()) || (*wh)->GetUncertainty().Z()==0){writableHits.erase(wh);}
       else++wh;
     }
   if(writableHits.size()>0){
   CP::THandle<CP::TReconCluster> wcl = CreateCluster("fortrack",writableHits.begin(),writableHits.end());
   
   writableClusters.push_back(wcl);
   }
  
 }
 writableHits.clear();
  }


	
  if(writableClusters.size()>0){
    
    CP::THandle<CP::TReconTrack> track(new CP::TReconTrack);
    track= CreateTrackFromClusters("TTracking3D_2",
				   writableClusters.begin(), writableClusters.end());
    if(track){
    match2.push_back(track);
    
    return true;}
  }
 
  return false;
  
}



void FindTrackCandidates(CP::TReconObjectContainer& tracksX,CP::TReconObjectContainer& tracksU,CP::TReconObjectContainer& tracksV,CP::TReconObjectContainer& match3,CP::TReconObjectContainer& match2){

    //********************************************************************************
  //Algorithm match 2D tracks by time(z) of their start and end positions
  //********************************************************************************

  double distCut=100;
 int trackNum=0;
  if(tracksX.size()>0 && tracksU.size()>0 && tracksV.size()>0) {

    
    
  CP::TReconObjectContainer::iterator trX = tracksX.begin();
  while(trX != tracksX.end()) {
    CP::THandle<CP::TReconTrack> trackX = *trX;

    CP::THandle<CP::THitSelection> hitX_t = trackX->GetHits();

    /* for(CP::THitSelection::iterator h = (*hitX_t).begin() ; h!=(*hitX_t).end();++h){
double ht=((*h)->GetConstituent()->GetTime()+1.6*unit::ms)/(500*unit::ns);
	    double hw=CP::GeomId::Captain::GetWireNumber((*h)->GetGeomId());
	    std::cout<<"xTIME="<<ht<<"; xWIRE="<<hw<<std::endl;
	    }*/


    
    double maxZX=MaxZ(trackX);
    double minZX=MinZ(trackX);
    int ntracks=0;
    double xuDiff=0;
    double xvDiff=0;
    double uvDiff=0;
    if(tracksU.size()>0){

      std::sort(tracksU.begin(),tracksU.end(),CompTracks(maxZX,minZX));
    ntracks++;
     xuDiff=abs(maxZX-MaxZ(tracksU[0]))+abs(minZX-MinZ(tracksU[0]));
    }
    if(tracksV.size()>0){

       std::sort(tracksV.begin(),tracksV.end(),CompTracks(maxZX,minZX));
    xvDiff=abs(maxZX-MaxZ(tracksV[0]))+abs(minZX-MinZ(tracksV[0]));
    ntracks++;
    }
   
    if(ntracks==2){
      //Turned ou that this dis should be 0 , otherwise nothing works :/ 
      double uvDiff=0;//abs(MaxZ(tracksV[0])-MaxZ(tracksU[0]))+abs(MinZ(tracksV[0])-MinZ(tracksU[0]));
        std::cout<<"uvDiff="<<uvDiff<<"; xuDiff="<<xuDiff<<"; xvDiff="<<xvDiff<<std::endl;
    if(uvDiff<distCut && xuDiff<distCut && xvDiff<distCut)
      {
	if(Assemble3DTrack(trackX,tracksU[0],tracksV[0],match3,trackNum)){
	  tracksX.erase(trX);
	  trackNum++;
	  tracksU.erase(tracksU.begin());
	  tracksV.erase(tracksV.begin());
	}else ++trX;
	}else ++trX;
    }
    if(ntracks==1 && xvDiff==0){
      std::cout<<"xuDiff="<<xuDiff<<std::endl;
       if(xuDiff<distCut )
      {
	if(Assemble2DTrack(trackX,tracksU[0],match2,trackNum)){ 
	tracksX.erase(trX);
	trackNum++;
	tracksU.erase(tracksU.begin());
	}else ++trX;
      }else ++trX;
    }
    if(ntracks==1 && xuDiff==0){
      std::cout<<"xvDiff="<<xvDiff<<std::endl;
       if(xvDiff<distCut )
      {
	if(Assemble2DTrack(trackX,tracksV[0],match2,trackNum)){ 
	tracksX.erase(trX);
	trackNum++;
	tracksU.erase(tracksU.begin());
	}else ++trX;
      }else ++trX;
    }
    if(ntracks==0)++trX;   
  }
  }
  if(tracksX.size()>0 && tracksU.size()>0) {
    for (CP::TReconObjectContainer::iterator trX = tracksX.begin();
       trX != tracksX.end(); ) {
    CP::THandle<CP::TReconTrack> trackX = *trX;
    double maxZX=MaxZ(trackX);
    double minZX=MinZ(trackX);
    double xuDiff=0;

    if(tracksU.size()>0){

       std::sort(tracksU.begin(),tracksU.end(),CompTracks(maxZX,minZX));
    xuDiff=abs(maxZX-MaxZ(tracksU[0]))+abs(minZX-MinZ(tracksU[0]));
  std::cout<<"xuDiff="<<xuDiff<<std::endl;
    if(xuDiff<distCut )
      {
	if(Assemble2DTrack(trackX,tracksU[0],match2,trackNum)){ 
	tracksX.erase(trX);
	trackNum++;
	tracksU.erase(tracksU.begin());
	}else ++trX;
      }else ++trX;
    }else ++trX; 
  }
  }
  if(tracksX.size()>0 && tracksV.size()>0) {
     for (CP::TReconObjectContainer::iterator trX = tracksX.begin();
       trX != tracksX.end(); ) {
    CP::THandle<CP::TReconTrack> trackX = *trX;
    double maxZX=MaxZ(trackX);
    double minZX=MinZ(trackX);
    double xvDiff=0;
    if(tracksV.size()>0){
       std::sort(tracksV.begin(),tracksV.end(),CompTracks(maxZX,minZX));
    double xvDiff=abs(maxZX-MaxZ(tracksV[0]))+abs(minZX-MinZ(tracksV[0]));
  std::cout<<"xuDiff="<<xvDiff<<std::endl;
    if(xvDiff<distCut)
      {
	if(Assemble2DTrack(trackX,tracksV[0],match2,trackNum)){
	tracksX.erase(trX);
	trackNum++;
	tracksV.erase(tracksV.begin());
	}else++trX;
      }else++trX;
    }else ++trX;
  }
  }


  

  
   if(tracksU.size()>0 && tracksV.size()>0) {
     for (CP::TReconObjectContainer::iterator trU = tracksU.begin();
       trU != tracksU.end(); ) {
    CP::THandle<CP::TReconTrack> trackU = *trU;
    double maxZU=MaxZ(trackU);
    double minZU=MinZ(trackU);
    double uvDiff=0;
    if(tracksV.size()>0){
       std::sort(tracksV.begin(),tracksV.end(),CompTracks(maxZU,minZU));
    double uvDiff=abs(maxZU-MaxZ(tracksV[0]))+abs(minZU-MinZ(tracksV[0]));
  std::cout<<"uvDiff="<<uvDiff<<std::endl;
    if(uvDiff<distCut)
      { 
	if(Assemble2DTrack(trackU,tracksV[0],match2,trackNum)){
	 
	tracksU.erase(trU);
	trackNum++;

	tracksV.erase(tracksV.begin());

	}else++trU;
      }else++trU;
    }else ++trU;
  }
  }
  
  
}
  

CP::TTracking3D::TTracking3D()
    : TAlgorithm("TTracking3D", 
                 "Break up objects into separate hits") {


 
}

CP::TTracking3D::~TTracking3D() { }



CP::THandle<CP::TAlgorithmResult>
CP::TTracking3D::Process(const CP::TAlgorithmResult& input,
                               const CP::TAlgorithmResult& input1,
                               const CP::TAlgorithmResult&) {

  //ReconObjects 
  CP::THandle<CP::TReconObjectContainer> inputObjects 
      = input.GetResultsContainer();
  
  //1
    CaptLog("TTracking3D Process " << GetEvent().GetContext());

    if (!inputObjects) {
        CaptError("No input objects");
        return CP::THandle<CP::TAlgorithmResult>();
    }
    
    // Create the output containers.
    CP::THandle<CP::TAlgorithmResult> result = CreateResult();
    std::unique_ptr<CP::TReconObjectContainer> 
        final(new CP::TReconObjectContainer("final"));


    //Copy inputObjects in "final" container
    CP::TReconObjectContainer tracksX;
    CP::TReconObjectContainer tracksU;
    CP::TReconObjectContainer tracksV;
    
    for (CP::TReconObjectContainer::iterator tr = inputObjects->begin();
         tr != inputObjects->end(); ++tr) {
        CP::THandle<CP::TReconTrack> track = *tr;
        if (!track) {
            final->push_back(*tr);
            continue;
        }
	int plane = CheckObjectPlaneT(track);
	if(plane==1){
	  tracksX.push_back(*tr);}
	else if(plane==2){
	  tracksU.push_back(*tr);}
	else if(plane==3){
	  tracksV.push_back(*tr);}
	else{std::cout<<"PLANEDEFININGFORCLUSTERSDOESNOTWORK"<<std::endl;}
    }

     CP::TReconObjectContainer match3;
    CP::TReconObjectContainer match2;
     CP::TReconObjectContainer used_clusters;
     FindTrackCandidates(tracksX,tracksU,tracksV,match3,match2);
     std::cout<<"MATCH3.size()="<<match3.size()<<std::endl;
     std::cout<<"MATCH2.size()="<<match2.size()<<std::endl;
     std::unique_ptr<CP::TReconObjectContainer> 
        match3Tr(new CP::TReconObjectContainer("match3Tr"));
     std::unique_ptr<CP::TReconObjectContainer> 
        match2Tr(new CP::TReconObjectContainer("match2Tr"));

     
     if(match3.size()>0)
       {
	 for(std::size_t i=0;i<match3.size();++i){
	   final->push_back(match3[i]);
	   match3Tr->push_back(match3[i]);
	 }
       }
     if(match2.size()>0)
       {
	 for(std::size_t i=0;i<match2.size();++i){
	   final->push_back(match2[i]);
	   match2Tr->push_back(match2[i]);
	 }
       }

       

    result->AddResultsContainer(match3Tr.release());
    result->AddResultsContainer(match2Tr.release());
    result->AddResultsContainer(final.release());

    return result;
}
