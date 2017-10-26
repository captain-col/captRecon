#ifndef TWireClustering_hxx_seen
#define TWireClustering_hxx_seen


#include <THandle.hxx>
#include <THitSelection.hxx>
#include <TChannelCalib.hxx>

typedef std::vector<CP::THandle<CP::TReconCluster>> clusterCollection;

namespace CP {
  class TWireClustering;
};

class CP::TWireClustering {
public:
  TWireClustering();
  clusterCollection WireCluster(int plane, CP::THitSelection hits,int minPts, double maxDistW, double maxDistT);
  
private:
  int fMinPoints;
  double fMaxDistW;
  double fMaxDistT;
  
};

CP::TWireClustering::TWireClustering(){};
clusterCollection CP::TWireClustering::WireCluster(int plane, CP::THitSelection hits,int minPts, double maxDistW, double maxDistT){

  typedef std::tuple<double,double,CP::THandle<CP::THit>> wTuple;
  
  clusterCollection wireClusters;

  double wireTimeStep = -1.0;
  double digitSampleOffset = 0;
  double digitSampleStep = -1;
  double triggerOffset = digitSampleOffset;
  CP::TChannelCalib chanCalib;
  
  std::vector<wTuple> wiresTimesCharges;
    
  for (CP::THitSelection::iterator h = hits.begin(); h != hits.end(); ++h) {


    //CP::THandle<CP::THit> test = (*h);
    
    CP::TGeometryId id = (*h)->GetGeomId();
    CP::TChannelId cid = (*h)->GetChannelId();
    if (CP::GeomId::Captain::GetWirePlane(id) != plane) continue;
    digitSampleStep = 0.5;
    if (wireTimeStep < 0) {
      wireTimeStep = chanCalib.GetTimeConstant(cid,1);
    }
    double timeUnit =  wireTimeStep/digitSampleStep;

    // The wire number (offset for the middle of the bin).
    double wire = CP::GeomId::Captain::GetWireNumber(id) + 0.5;
    // The hit charge
    //double charge = (*h)->GetCharge();
    // The hit time.
    double hTime = (*h)->GetTime();
    hTime = hTime;
    // The digitized hit time.
    double dTime = hTime/timeUnit + triggerOffset;

    //std::cout<<"Time="<<dTime<<std::endl;
    
    wiresTimesCharges.push_back(std::make_tuple(wire,dTime,(*h)));	            
  }

  std::sort(wiresTimesCharges.begin(), wiresTimesCharges.end(), [] (wTuple const& a, wTuple const& b) { return std::get<0>(a) < std::get<0>(b); });


  std::vector<std::vector<wTuple>> trackSeeds;
  std::vector<std::vector<wTuple>> tracks;
  std::vector<double> hotWires;
  double hotWire = 0;
  int nhWires = 0;
    
  for (auto w:wiresTimesCharges) {
    std::vector<wTuple> track_tmp;
    track_tmp.push_back(w);
    trackSeeds.push_back(track_tmp);
    if(std::get<0>(w) == hotWire) {
      nhWires++;
    }
    else{
      if (nhWires > 20) hotWires.push_back(hotWire);
      nhWires = 0;
    }     
    if (w == wiresTimesCharges.back()) {
      if (nhWires > 20) hotWires.push_back(hotWire);
    }
    hotWire = std::get<0>(w);
  }
  int nTracks = 0;
  //std::cout<<"Plane "<<planes[plane]<<std::endl;
  std::cout<<"wires="<<wiresTimesCharges.size()<<std::endl;   
  // Veto hot wires
  std::vector<int> hitsToVeto;
  for (auto hw:hotWires) {
    for (unsigned int i = 0; i < trackSeeds.size(); i++) {
      if (hw == std::get<0>(trackSeeds[i][0])) {
	hitsToVeto.push_back(i);
      }
    }
  }

  reverse(hitsToVeto.begin(),hitsToVeto.end());
    
  for (auto i:hitsToVeto) {
    trackSeeds.erase(trackSeeds.begin() + i);
    wiresTimesCharges.erase(wiresTimesCharges.begin() + i);
  }
  std::cout<<"track seeds="<<trackSeeds.size()<<std::endl;
        
  for (unsigned int i = 0; i < trackSeeds.size(); i++) {
    for (unsigned int k = 0; k < trackSeeds[i].size(); k++) {	
      for (int j = wiresTimesCharges.size() -1; j  >= 0; j--) {	  
	double deltaW = fabs(std::get<0>(trackSeeds[i][k])-std::get<0>(wiresTimesCharges[j]));
	double deltaT = fabs(std::get<1>(trackSeeds[i][k])-std::get<1>(wiresTimesCharges[j]));
	//std::cout<<"Deltas="<<deltaW<<" "<<deltaT<<std::endl;
	if ((deltaW !=0 && deltaT != 0) && deltaW < maxDistW && deltaT < maxDistT) {
	  trackSeeds[i].push_back(wiresTimesCharges[j]);
	  wiresTimesCharges.erase(wiresTimesCharges.begin() + j);
	}	  
      }	
    }
  }
  std::cout<<"pruned seeds="<<trackSeeds.size()<<std::endl;
  for (auto ts:trackSeeds) {
    if (ts.size() > minPts) {
      tracks.push_back(ts);
      nTracks++;
    }
  }
  
  std::cout<<"ntracks = "<<nTracks<<std::endl;
  for (auto t:tracks) {
    CP::THitSelection hitContainer;
    for (auto ts:t) {
      hitContainer.AddHit(std::get<2>(ts));
    }
    
    CP::THandle<CP::TReconCluster> tmpCluster
      = CreateCluster("tmpCluster",hitContainer.begin(),hitContainer.end());      

    wireClusters.push_back(tmpCluster);
    
  }
  
  
  return wireClusters;

}


#endif
