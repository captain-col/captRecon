#include "TTracking3D.hxx"
#include "HitUtilities.hxx"
#include "CompareReconObjects.hxx"
#include "CreateTrack.hxx"
#include "CreateHit.hxx"

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

  int CheckObjectPlaneT(const CP::THandle<CP::TReconBase>& obj){
    
 
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

bool CompHitsZ(const CP::THandle<CP::THit>& lhs, const CP::THandle<CP::THit>& rhs){
  return lhs->GetPosition().Z() < rhs->GetPosition().Z();
}

bool CheckPosition(TVector3 s, TVector3 f, TVector3 posX){
  if(posX.Z()>s.Z() && posX.Z()<f.Z()) return true;
   if(posX.Z()<s.Z() && posX.Z()>f.Z()) return true;
  return false;
}

TVector3 Intersection(TVector3 pos1S,TVector3 pos1F,TVector3 pos2S,TVector3 pos2F){

  TVector3 intersect;

  if(pos1S.Y()==pos1F.Y() ){
    double b = std::rand()% 3+1;
    b/=5000;
    pos1S.SetY(b);
  }
  if(pos2S.Y()==pos2F.Y() ){
    double b = std::rand()% 3+1;
    b/=5000;
    pos2S.SetY(b);
  }
  //define initial lines

  double a1=0;
  double a2=0;

  a1=(pos1S.Y()-pos1F.Y())/(pos1S.X()-pos1F.X());
  a2=(pos2S.Y()-pos2F.Y())/(pos2S.X()-pos2F.X());
  std::cout<<"a1="<<a1<<"; a2="<<a2<<std::endl;

  //define ortogonal lines through 1st given point

  double m1=0;
  double c1=0;
  double m2=0;
  double c2=0;

  m1=-1/a1;
  c1=pos1S.Y()+m1*pos1S.X();
  m2=-1/a2;
  c2=pos2S.Y()+m2*pos2S.X();
  std::cout<<"m1="<<m1<<"; m2="<<m2<<std::endl;
    std::cout<<"c1="<<c1<<"; c2="<<c2<<std::endl;
  //find intersection

  intersect.SetX((c2-c1)/(m1-m2));
  intersect.SetY(m1*(c2-c1)/(m1-m2)+c1);
  intersect.SetZ(0.0);
  

  


  return intersect;
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

bool Assemble3DTrack( CP::THandle<CP::TReconTrack> trackX, CP::THandle<CP::TReconTrack> trackU, CP::THandle<CP::TReconTrack> trackV,CP::TReconObjectContainer& match3){
  CP::THandle<CP::THitSelection> hitX = trackX->GetHits();
  CP::THandle<CP::THitSelection> hitU = trackU->GetHits();
  CP::THandle<CP::THitSelection> hitV = trackV->GetHits();
  
  //Sort Hits by wire number min to max
  std::sort((*hitX).begin(),(*hitX).end(),CompHitsZ);
  std::sort((*hitU).begin(),(*hitU).end(),CompHitsZ);
  std::sort((*hitV).begin(),(*hitV).end(),CompHitsZ);
  
  TVector3 frontX = (*hitX).front()->GetPosition();
  TVector3 frontU = (*hitU).front()->GetPosition();
  TVector3 frontV = (*hitV).front()->GetPosition();
  TVector3 backX = (*hitX).back()->GetPosition();
  TVector3 backU = (*hitU).back()->GetPosition();
  TVector3 backV = (*hitV).back()->GetPosition();

 TVector3 interXU_begin=Intersection(frontX,backX,frontU,backU);
  TVector3 interXU_end=Intersection(backX,frontX,backU,frontU);

  TVector3 interXV_begin=Intersection(frontX,backX,frontV,backV);
  TVector3 interXV_end=Intersection(backX,frontX,backV,frontV);

  TVector3 startPoint;
  TVector3 endPoint;
 
  if((interXU_begin-interXV_begin).Mag()<100 && (interXU_end-interXV_end).Mag()<100){
    startPoint.SetX((interXU_begin.X()+interXV_begin.X())/2);
    startPoint.SetY((interXU_begin.Y()+interXV_begin.Y())/2);
    startPoint.SetZ(frontX.Z());
    endPoint.SetX((interXU_end.X()+interXV_end.X())/2);
    endPoint.SetY((interXU_end.Y()+interXV_end.Y())/2);
    endPoint.SetZ(backX.Z());
  }
  else
    {
      startPoint.SetX(interXU_begin.X());
    startPoint.SetY(interXU_begin.Y());
    startPoint.SetZ(frontX.Z());
    endPoint.SetX(interXU_end.X());
    endPoint.SetY(interXU_end.Y());
    endPoint.SetZ(backX.Z()); 
    }
  
  double tLength = (startPoint-endPoint).Mag();

  int nparts = (int)tLength/3;

  if(nparts<1) nparts=1;
  if(nparts > 20 && nparts<51) nparts/=2;
  if(nparts > 50 && nparts<101) nparts/=5;
  if(nparts > 100) nparts/=10;
  TVector3 stepVect((endPoint.X()-startPoint.X())/nparts,(endPoint.Y()-startPoint.Y())/nparts,(endPoint.Z()-startPoint.Z())/nparts);

      CP::TReconObjectContainer writableClusters;
  for(std::size_t l=0;l<nparts;++l){

    
    TVector3 s = startPoint + l*stepVect;
    TVector3 f = startPoint + (l+1)*stepVect;

    CP::THitSelection writableHits;
 CP::THitSelection inHitsX;
 CP::THitSelection inHitsU;
 CP::THitSelection inHitsV;
 for(CP::THitSelection::iterator hx = (*hitX).begin(); hx!=(*hitX).end();++hx){
   CP::THandle<CP::THit> hittX = (*hx);
   TVector3 posX = hittX->GetPosition();

   if(CheckPosition(s,f,posX)){
     inHitsX.push_back(hittX);
   }
 }
 
 for(CP::THitSelection::iterator hu = (*hitU).begin(); hu!=(*hitU).end();++hu){
   CP::THandle<CP::THit> hittU = (*hu);
   TVector3 posU = hittU->GetPosition();
   if(CheckPosition(s,f,posU)){
     inHitsU.push_back(hittU);
   }
 }
 for(CP::THitSelection::iterator hv = (*hitV).begin(); hv!=(*hitV).end();++hv){
   CP::THandle<CP::THit> hittV = (*hv);
   TVector3 posV = hittV->GetPosition();
   if(CheckPosition(s,f,posV)){
     inHitsV.push_back(hittV);
   }
 }
 
 if(inHitsX.size()>0){
   
   for(std::size_t i=0;i<inHitsX.size();++i){
 
   if(inHitsU.size()>0){
     for(std::size_t j=0;j<inHitsU.size();++j){

       if(inHitsV.size()>0){
	 for(std::size_t k=0;k<inHitsV.size();++k){

	   CP::THandle<CP::THit> h1=inHitsX[i]->GetConstituent();
	   CP::THandle<CP::THit> h2=inHitsU[j]->GetConstituent();
	   CP::THandle<CP::THit> h3=inHitsV[k]->GetConstituent();
	   TVector3 hitPosition;
	 
	   OverlapedXY(h1,h2,h3,hitPosition);
	  
	   if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
	 }
       }
       CP::THandle<CP::THit> h1=inHitsX[i]->GetConstituent();
       CP::THandle<CP::THit> h2=inHitsU[j]->GetConstituent();
       if(!CreateHit(writableHits,PositionXY(h1,h2),0.0,h1,h2,CP::THandle<CP::THit>()))continue;
     }
   }
   if(inHitsV.size()>0){
     for(std::size_t k=0;k<inHitsV.size();++k){
   
       CP::THandle<CP::THit> h1=inHitsX[i]->GetConstituent();
       CP::THandle<CP::THit> h2=inHitsV[k]->GetConstituent();
       if(!CreateHit(writableHits,PositionXY(h1,h2),0.0,h1,h2,CP::THandle<CP::THit>()))continue;
     }
   }
   }
 }
 
 if(writableHits.size()>0){
   CP::THandle<CP::TReconCluster> wcl = CreateCluster("fortrack",writableHits.begin(),writableHits.end());
   writableClusters.push_back(wcl);
 writableHits.clear();
  }
  }
 
  if(writableClusters.size()>0){
    CP::THandle<CP::TReconTrack> track(new CP::TReconTrack);
    track= CreateTrackFromClusters("TTracking3D",
				   writableClusters.begin(), writableClusters.end());

    match3.push_back(track);

    return true;
  }
  return false;
  
}

bool Assemble2DTrack( CP::THandle<CP::TReconTrack> trackX, CP::THandle<CP::TReconTrack> trackU,CP::TReconObjectContainer& match2){

  CP::THandle<CP::THitSelection> hitX = trackX->GetHits();
  CP::THandle<CP::THitSelection> hitU = trackU->GetHits();

  //Sort Hits by wire number min to max
  std::sort((*hitX).begin(),(*hitX).end(),CompHitsZ);
  std::sort((*hitU).begin(),(*hitU).end(),CompHitsZ);


   TVector3 frontX = (*hitX).front()->GetPosition();
  TVector3 frontU = (*hitU).front()->GetPosition();

  TVector3 backX = (*hitX).back()->GetPosition();
  TVector3 backU = (*hitU).back()->GetPosition();


   if(((frontX.X()-backX.X())==0) || ((frontX.Y()-backX.Y())==0) || ((frontU.X()-backU.X())==0) || ((frontU.Y()-backU.Y())==0)) return false;

 TVector3 interXU_begin=Intersection(frontX,backX,frontU,backU);
  TVector3 interXU_end=Intersection(backX,frontX,backU,frontU);

  TVector3 startPoint;
  TVector3 endPoint;
  
    startPoint.SetX(interXU_begin.X());
    startPoint.SetY(interXU_begin.Y());
    startPoint.SetZ((*(*hitX).begin())->GetPosition().Z());
    endPoint.SetX(interXU_end.X());
    endPoint.SetY(interXU_end.Y());
    endPoint.SetZ((*(*hitX).end())->GetPosition().Z());

  
  double tLength = (startPoint-endPoint).Mag();
  int nparts = (int)tLength/3;
  if(nparts > 20) nparts/=2;
  if(nparts > 50) nparts/=5;
  if(nparts > 100) nparts/=10;
  TVector3 stepVect((endPoint.X()-startPoint.X())/nparts,(endPoint.Y()-startPoint.Y())/nparts,(endPoint.Z()-startPoint.Z())/nparts);

  std::vector<CP::THandle<CP::THitSelection>> parts;
      CP::TReconObjectContainer writableClusters;
  for(std::size_t i=0;i<nparts;++i){

    TVector3 s = startPoint + i*stepVect;
    TVector3 f = startPoint + (i+1)*stepVect;
 CP::THitSelection writableHits;
 CP::THitSelection inHitsX;
 CP::THitSelection inHitsU;
 CP::THitSelection inHitsV;
 for(CP::THitSelection::iterator hx = (*hitX).begin(); hx!=(*hitX).end();++hx){
   CP::THandle<CP::THit> hittX = (*hx);
   TVector3 posX = hittX->GetPosition();
   if(CheckPosition(s,f,posX)){
     inHitsX.push_back(hittX);
   }
 }
 
 for(CP::THitSelection::iterator hu = (*hitU).begin(); hu!=(*hitU).end();++hu){
   CP::THandle<CP::THit> hittU = (*hu);
   TVector3 posU = hittU->GetPosition();
   if(CheckPosition(s,f,posU)){
     inHitsU.push_back(hittU);
   }
 }


 if(inHitsX.size()>0){
   TVector3 hitPosition;
   for(std::size_t i=0;i<inHitsX.size();++i){
     
     if(inHitsU.size()>0){
     for(std::size_t k=0;k<inHitsU.size();++k){
       CP::THandle<CP::THit> h1=inHitsX[i]->GetConstituent();
       CP::THandle<CP::THit> h2=inHitsU[k]->GetConstituent();
       if(!CreateHit(writableHits,PositionXY(h1,h2),0.0,h1,h2,CP::THandle<CP::THit>()))continue;
     }
   }
   }
 }
 
 if(writableHits.size()>0){
   CP::THandle<CP::TReconCluster> wcl = CreateCluster("fortrack",writableHits.begin(),writableHits.end());
   writableClusters.push_back(wcl);}
 writableHits.clear();
  }

  if(writableClusters.size()>0){
    CP::THandle<CP::TReconTrack> track(new CP::TReconTrack);
    track= CreateTrackFromClusters("TTracking3D",
				   writableClusters.begin(), writableClusters.end());
    match2.push_back(track);
    return true;
  }
  
  return false;
  
}



void FindTrackCandidates(CP::TReconObjectContainer& tracksX,CP::TReconObjectContainer& tracksU,CP::TReconObjectContainer& tracksV,CP::TReconObjectContainer& match3,CP::TReconObjectContainer& match2){

  if(tracksX.size()>0 && tracksU.size()>0 && tracksV.size()>0) {
  CP::TReconObjectContainer::iterator trX = tracksX.begin();
  while(trX != tracksX.end()) {
    CP::THandle<CP::TReconTrack> trackX = *trX;
    double maxZX=MaxZ(trackX);
    double minZX=MinZ(trackX);

    std::sort(tracksU.begin(),tracksU.end(),[maxZX,minZX](const CP::THandle<CP::TReconTrack>& lhs, const CP::THandle<CP::TReconTrack>& rhs)
	      {
		double l=abs(MaxZ(lhs)-maxZX)+abs(MinZ(lhs)-minZX);
		double r=abs(MaxZ(rhs)-maxZX)+abs(MinZ(rhs)-minZX);
		return l<r;});
    std::sort(tracksV.begin(),tracksV.end(),[maxZX,minZX](const CP::THandle<CP::TReconTrack>& lhs, const CP::THandle<CP::TReconTrack>& rhs)
	      {
		double l=abs(MaxZ(lhs)-maxZX)+abs(MinZ(lhs)-minZX);
		double r=abs(MaxZ(rhs)-maxZX)+abs(MinZ(rhs)-minZX);
		return l<r;});
    double uvDiff=abs(MaxZ(tracksV[0])-MaxZ(tracksU[0]))+abs(MinZ(tracksV[0])-MinZ(tracksU[0]));
    double xuDiff=abs(maxZX-MaxZ(tracksU[0]))+abs(minZX-MinZ(tracksU[0]));
    double xvDiff=abs(maxZX-MaxZ(tracksV[0]))+abs(minZX-MinZ(tracksV[0]));
    std::cout<<"uvDiff="<<uvDiff<<"; xuDiff="<<xuDiff<<"; xvDiff="<<xvDiff<<std::endl;
    if(uvDiff<100 && xuDiff<100 && xvDiff<100)
      {
	if(Assemble3DTrack(trackX,tracksU[0],tracksV[0],match3)){
	  tracksX.erase(trX);
	  std::cout<<"ERASE"<<std::endl;
	  tracksU.erase(tracksU.begin());
	  tracksV.erase(tracksV.begin());
	}else ++trX;
	}else ++trX;

  }
  }
  if(tracksX.size()>0 && tracksU.size()>0 && tracksV.size()==0) {
    for (CP::TReconObjectContainer::iterator trX = tracksX.begin();
       trX != tracksX.end(); ) {
    CP::THandle<CP::TReconTrack> trackX = *trX;
    double maxZX=MaxZ(trackX);
    double minZX=MinZ(trackX);
    
    std::sort(tracksU.begin(),tracksU.end(),[maxZX,minZX](const CP::THandle<CP::TReconTrack>& lhs, const CP::THandle<CP::TReconTrack>& rhs)
	      {
		double l=abs(MaxZ(lhs)-maxZX)+abs(MinZ(lhs)-minZX);
		double r=abs(MaxZ(rhs)-maxZX)+abs(MinZ(rhs)-minZX);
		return l<r;});
    double xuDiff=abs(maxZX-MaxZ(tracksU[0]))+abs(minZX-MinZ(tracksU[0]));
 
    if(xuDiff<70 )
      {
	if(Assemble2DTrack(trackX,tracksU[0],match2)){ 
	tracksX.erase(trX);
	tracksU.erase(tracksU.begin());
	}else ++trX;
      }else ++trX;
      
  }
  }
  if(tracksX.size()>0 && tracksV.size()>0 && tracksU.size()==0) {
     for (CP::TReconObjectContainer::iterator trX = tracksX.begin();
       trX != tracksX.end(); ) {
    CP::THandle<CP::TReconTrack> trackX = *trX;
    double maxZX=MaxZ(trackX);
    double minZX=MinZ(trackX);
    
    std::sort(tracksV.begin(),tracksV.end(),[maxZX,minZX](const CP::THandle<CP::TReconTrack>& lhs, const CP::THandle<CP::TReconTrack>& rhs)
	      {
		double l=abs(MaxZ(lhs)-maxZX)+abs(MinZ(lhs)-minZX);
		double r=abs(MaxZ(rhs)-maxZX)+abs(MinZ(rhs)-minZX);
		return l<r;});
    double xvDiff=abs(maxZX-MaxZ(tracksV[0]))+abs(minZX-MinZ(tracksV[0]));
 
    if(xvDiff<70)
      {
	if(Assemble2DTrack(trackX,tracksV[0],match2)){
	tracksX.erase(trX);
	tracksV.erase(tracksV.begin());
	}else++trX;
      }else++trX;
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
