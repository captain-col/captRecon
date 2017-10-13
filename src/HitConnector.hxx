
#ifndef HitConnector_hxx_seen
#define HitConnector_hxx_seen

#include "HitUtilities.hxx"
#include "ECaptRecon.hxx"

#include <TDriftPosition.hxx>
#include <TUnitsTable.hxx>
#include <THandle.hxx>
#include <TReconHit.hxx>
#include <THit.hxx>
#include <CaptGeomId.hxx>
#include <TRuntimeParameters.hxx>


/*namespace CP{

  CP::THitSelection HitConnector3D(const CP::THandle<CP::THitSelection>& xHits, const CP::THandle<CP::THitSelection>& uHits, const CP::THandle<CP::THitSelection>& vHits);

  CP::THitSelection HitConnector2D(const CP::THandle<CP::THitSelection>& xHits, const CP::THandle<CP::THitSelection>& uHits);

  };*/


// This algoright takes 3(2) arrays of already sorted 2D Hits and form array of //3D hits such, that only nearby hits in each arrays are connected. For example //0th heat in array 1 is connected with 0 and 1st hits in array two. Hit ith is //connected with i-1,i,i+1 hits in over array.

bool CompHitsX(const CP::THandle<CP::THit>& lhs, const CP::THandle<CP::THit>& rhs){
  return lhs->GetPosition().X() < rhs->GetPosition().X();
}


CP::THitSelection HitConnector3D(const CP::THandle<CP::THitSelection>& xHits, const CP::THandle<CP::THitSelection>& uHits, const CP::THandle<CP::THitSelection>& vHits){

  int sizeX = xHits->size();
  int sizeU = uHits->size();
  int sizeV = vHits->size();
  int maxS=std::max(sizeX,sizeU);
  maxS=std::max(maxS,sizeV);
  int minS=std::min(sizeX,sizeU);
  minS=std::min(minS,sizeV);
  int order1 = 1;
  int order2 = 2;
  int order3 = 3;
  CP::THitSelection hits1;
  CP::THitSelection hits2;
  CP::THitSelection hits3;
  if(maxS==sizeX)order1=1;
  if(maxS==sizeU)order1=2;
  if(maxS==sizeV)order1=3;
  if(minS==sizeX)order3=1;
  if(minS==sizeU)order3=2;
  if(minS==sizeV)order3=3;
  if(order1==1 && order3==2)order2=3;
  if(order1==2 && order3==1)order2=3;
  if(order1==1 && order3==3)order2=2;
  if(order1==3 && order3==1)order2=2;
  if(order1==2 && order3==3)order2=1;
  if(order1==3 && order3==2)order2=1;
  int orderF=1;
  if(order1==1 && order2==2 && order3==3){
    hits1=(*xHits);
    hits2=(*uHits);
    hits3=(*vHits);
	    }

  if(order1==1 && order2==3 && order3==2){
    hits1=(*xHits);
    hits2=(*vHits);
    hits3=(*uHits);
	    }
  if(order1==2 && order2==1 && order3==3){
    hits1=(*uHits);
    hits2=(*xHits);
    hits3=(*vHits);
    orderF=2;
	    }
  if(order1==2 && order2==3 && order3==1){
    hits1=(*uHits);
    hits2=(*vHits);
    hits3=(*xHits);
    orderF=3;
	    }
  if(order1==3 && order2==2 && order3==1){
    hits1=(*vHits);
    hits2=(*uHits);
    hits3=(*xHits);
    orderF=3;
	    }
  if(order1==3 && order2==1 && order3==2){
    hits1=(*vHits);
    hits2=(*xHits);
    hits3=(*uHits);
    orderF=2;
	    }
  CP::THitSelection writableHits;
  // std::cout<<"1Size="<<hits1.size()<<"; 2Size="<<hits2.size()<<"; 3Size="<<hits3.size()<<std::endl;
  for(int i=0;i<hits1.size();++i){
    int dif_i = hits2.size()-i;
    if(i==0){
      if(dif_i>1){
	for(int j=0;j<2;++j){
	  int dif_j = hits3.size()-j;
	  if(j==0){
	    if(dif_j>1){
	      for(int k=0;k<2;++k){
		if(orderF==1){
		  CP::THandle<CP::THit> h1=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==2){
		  CP::THandle<CP::THit> h2=hits1[i];
		  CP::THandle<CP::THit> h1=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==3){
		  CP::THandle<CP::THit> h3=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h1=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
	      }
	    }else{
	      for(int k=0;k<hits3.size();++k){
		if(orderF==1){
		  CP::THandle<CP::THit> h1=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==2){
		  CP::THandle<CP::THit> h2=hits1[i];
		  CP::THandle<CP::THit> h1=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==3){
		  CP::THandle<CP::THit> h3=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h1=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
	      }
	    }
	  }else{
	    if(dif_j>2){
	      for(int k=j-1;k<j+2;++k){
		if(orderF==1){
		  CP::THandle<CP::THit> h1=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==2){
		  CP::THandle<CP::THit> h2=hits1[i];
		  CP::THandle<CP::THit> h1=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==3){
		  CP::THandle<CP::THit> h3=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h1=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
	      }
	    }else{
	      int bound = 0;
	      if(hits3.size()==1)bound=1;
	      if(hits3.size()==2)bound=2;
	      if(hits3.size()>2)bound=3;
	      bound = hits3.size()-bound;
	      for(int k=bound;k<hits3.size();++k){
		if(orderF==1){
		  CP::THandle<CP::THit> h1=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==2){
		  CP::THandle<CP::THit> h2=hits1[i];
		  CP::THandle<CP::THit> h1=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==3){
		  CP::THandle<CP::THit> h3=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h1=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
	      }
	    }
	  }
	}
      }else{//hits2.size-i <=1		
	for(int j=0;j<hits2.size();++j){
	  int dif_j = hits3.size()-j;
	  if(j==0){
	    if(dif_j>1){
	      for(int k=0;k<2;++k){
		if(orderF==1){
		  CP::THandle<CP::THit> h1=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==2){
		  CP::THandle<CP::THit> h2=hits1[i];
		  CP::THandle<CP::THit> h1=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==3){
		  CP::THandle<CP::THit> h3=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h1=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
	      }
	    }else{
	      for(int k=0;k<hits3.size();++k){
		if(orderF==1){
		  CP::THandle<CP::THit> h1=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==2){
		  CP::THandle<CP::THit> h2=hits1[i];
		  CP::THandle<CP::THit> h1=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==3){
		  CP::THandle<CP::THit> h3=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h1=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
	      }
	    }
	  }else{
	    if(dif_j>2){
	      for(int k=j-1;k<j+2;++k){
		if(orderF==1){
		  CP::THandle<CP::THit> h1=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==2){
		  CP::THandle<CP::THit> h2=hits1[i];
		  CP::THandle<CP::THit> h1=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==3){
		  CP::THandle<CP::THit> h3=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h1=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
	      }
	    }else{
	      int bound = 0;
	      if(hits3.size()==1)bound=1;
	      if(hits3.size()==2)bound=2;
	      if(hits3.size()>2)bound=3;
	      bound = hits3.size()-bound;
	      for(int k=bound;k<hits3.size();++k){
		if(orderF==1){
		  CP::THandle<CP::THit> h1=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==2){
		  CP::THandle<CP::THit> h2=hits1[i];
		  CP::THandle<CP::THit> h1=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==3){
		  CP::THandle<CP::THit> h3=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h1=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
	      }
	    }
	  }
	}	
      }
    }else{  // if i > 0
 if(dif_i>2){
	for(int j=i-1;j<i+2;++j){
	  int dif_j = hits3.size()-j;
	  if(j==0){
	    if(dif_j>1){
	      for(int k=0;k<2;++k){
		if(orderF==1){
		  CP::THandle<CP::THit> h1=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==2){
		  CP::THandle<CP::THit> h2=hits1[i];
		  CP::THandle<CP::THit> h1=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==3){
		  CP::THandle<CP::THit> h3=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h1=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
	      }
	    }else{
	      for(int k=0;k<hits3.size();++k){
		if(orderF==1){
		  CP::THandle<CP::THit> h1=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==2){
		  CP::THandle<CP::THit> h2=hits1[i];
		  CP::THandle<CP::THit> h1=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==3){
		  CP::THandle<CP::THit> h3=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h1=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
	      }
	    }
	  }else{
	    if(dif_j>2){
	      for(int k=j-1;k<j+2;++k){
		if(orderF==1){
		  CP::THandle<CP::THit> h1=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==2){
		  CP::THandle<CP::THit> h2=hits1[i];
		  CP::THandle<CP::THit> h1=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==3){
		  CP::THandle<CP::THit> h3=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h1=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
	      }
	    }else{
	      int bound = 0;
	      if(hits3.size()==1)bound=1;
	      if(hits3.size()==2)bound=2;
	      if(hits3.size()>2)bound=3;
	      bound = hits3.size()-bound;
	      for(int k=bound;k<hits3.size();++k){
		if(orderF==1){
		  CP::THandle<CP::THit> h1=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==2){
		  CP::THandle<CP::THit> h2=hits1[i];
		  CP::THandle<CP::THit> h1=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==3){
		  CP::THandle<CP::THit> h3=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h1=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
	      }
	    }
	  }
	}
 }else{//hits2.size-i <=1
   int boundj = 0;
	      if(hits2.size()==1)boundj=1;
	      if(hits2.size()==2)boundj=2;
	      if(hits2.size()>2)boundj=3;
	      boundj = hits2.size()-boundj;
	      for(int j=boundj;j<hits2.size();++j){
		      int dif_j = hits3.size()-j;
	  if(j==0){
	    if(dif_j>1){
	      for(int k=0;k<2;++k){
		if(orderF==1){
		  CP::THandle<CP::THit> h1=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==2){
		  CP::THandle<CP::THit> h2=hits1[i];
		  CP::THandle<CP::THit> h1=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==3){
		  CP::THandle<CP::THit> h3=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h1=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
	      }
	    }else{
	      for(int k=0;k<hits3.size();++k){
		if(orderF==1){
		  CP::THandle<CP::THit> h1=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==2){
		  CP::THandle<CP::THit> h2=hits1[i];
		  CP::THandle<CP::THit> h1=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==3){
		  CP::THandle<CP::THit> h3=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h1=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
	      }
	    }
	  }else{
	    if(dif_j>2){
	      for(int k=j-1;k<j+2;++k){
		if(orderF==1){
		  CP::THandle<CP::THit> h1=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==2){
		  CP::THandle<CP::THit> h2=hits1[i];
		  CP::THandle<CP::THit> h1=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==3){
		  CP::THandle<CP::THit> h3=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h1=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
	      }
	    }else{
	      int bound = 0;
	      if(hits3.size()==1)bound=1;
	      if(hits3.size()==2)bound=2;
	      if(hits3.size()>2)bound=3;
	      bound = hits3.size()-bound;
	      for(int k=bound;k<hits3.size();++k){
		if(orderF==1){
		  CP::THandle<CP::THit> h1=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==2){
		  CP::THandle<CP::THit> h2=hits1[i];
		  CP::THandle<CP::THit> h1=hits2[j];
		  CP::THandle<CP::THit> h3=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
		if(orderF==3){
		  CP::THandle<CP::THit> h3=hits1[i];
		  CP::THandle<CP::THit> h2=hits2[j];
		  CP::THandle<CP::THit> h1=hits3[k];
		  TVector3 hitPosition;
		  OverlapedXY(h1,h2,h3,hitPosition);
		  if(!CreateHit(writableHits,hitPosition,0.0,h1,h2,h3))continue;
		}
	      }
	    }
	  }
	}	
 }
    }//else for i>0

  } //i cycle
    std::cout<<"#3DHitsCreated"<<writableHits.size()<<std::endl;
    //   std::sort(writableHits.begin(),writableHits.end(),CompHitsX);
  return writableHits;
}












CP::THitSelection HitConnector2D(const CP::THandle<CP::THitSelection>& xHits, const CP::THandle<CP::THitSelection>& uHits){

  int sizeX = xHits->size();
  int sizeU = uHits->size();

  int maxS=std::max(sizeX,sizeU);


  CP::THitSelection hits1;
  CP::THitSelection hits2;

  int orderF=1;
  if(maxS==sizeX){
    hits1=(*xHits);
    hits2=(*uHits);
	    }

  if(maxS==sizeU){
    hits1=(*uHits);
    hits2=(*xHits);
    orderF=2;
	    }


CP::THitSelection writableHits;
  for(int i=0;i<hits1.size();++i){
    int dif_i = hits2.size()-i;
    if(i==0){
      if(dif_i>1){
	for(int j=0;j<2;++j){
	  if(orderF==1){
	    CP::THandle<CP::THit> h1=hits1[i];
	    CP::THandle<CP::THit> h2=hits2[j];
	    if(!CreateHit(writableHits,PositionXY(h1,h2),0.0,h1,h2,CP::THandle<CP::THit>()))continue;
	  }
	  if(orderF==2){
	    CP::THandle<CP::THit> h2=hits1[i];
	    CP::THandle<CP::THit> h1=hits2[j];
	    if(!CreateHit(writableHits,PositionXY(h1,h2),0.0,h1,h2,CP::THandle<CP::THit>()))continue;
	  }
	}
      }else{//hits2.size-i <=1 
	for(int j=0;j<hits2.size();++j){
	  	  if(orderF==1){
	    CP::THandle<CP::THit> h1=hits1[i];
	    CP::THandle<CP::THit> h2=hits2[j];
	    if(!CreateHit(writableHits,PositionXY(h1,h2),0.0,h1,h2,CP::THandle<CP::THit>()))continue;
	  }
	  if(orderF==2){
	    CP::THandle<CP::THit> h2=hits1[i];
	    CP::THandle<CP::THit> h1=hits2[j];
	    if(!CreateHit(writableHits,PositionXY(h1,h2),0.0,h1,h2,CP::THandle<CP::THit>()))continue;
	  }

	}	
      }
    }else{  // if i > 0


 if(dif_i>2){
	for(int j=i-1;j<i+2;++j){
	  	  if(orderF==1){
	    CP::THandle<CP::THit> h1=hits1[i];
	    CP::THandle<CP::THit> h2=hits2[j];
	    if(!CreateHit(writableHits,PositionXY(h1,h2),0.0,h1,h2,CP::THandle<CP::THit>()))continue;
	  }
	  if(orderF==2){
	    CP::THandle<CP::THit> h2=hits1[i];
	    CP::THandle<CP::THit> h1=hits2[j];
	    if(!CreateHit(writableHits,PositionXY(h1,h2),0.0,h1,h2,CP::THandle<CP::THit>()))continue;
	  }

	}
 }else{//hits2.size-i <=1
   int boundj = 0;
	      if(hits2.size()==1)boundj=1;
	      if(hits2.size()==2)boundj=2;
	      if(hits2.size()>2)boundj=3;
	      boundj = hits2.size()-boundj;
	      for(int j=boundj;j<hits2.size();++j){

			  if(orderF==1){
	    CP::THandle<CP::THit> h1=hits1[i];
	    CP::THandle<CP::THit> h2=hits2[j];
	    if(!CreateHit(writableHits,PositionXY(h1,h2),0.0,h1,h2,CP::THandle<CP::THit>()))continue;
	  }
	  if(orderF==2){
	    CP::THandle<CP::THit> h2=hits1[i];
	    CP::THandle<CP::THit> h1=hits2[j];
	    if(!CreateHit(writableHits,PositionXY(h1,h2),0.0,h1,h2,CP::THandle<CP::THit>()))continue;
	  }

	}	
 }
    }//else for i>0

  }//i cycle
  // std::sort(writableHits.begin(),writableHits.end(),CompHitsX);
   return writableHits;

}



#endif


