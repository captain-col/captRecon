#ifndef CreateHit_hxx_seen
#define CreateHit_hxx_seen

#include "HitUtilities.hxx"
#include "ECaptRecon.hxx"

#include <TReconHit.hxx>
#include <TDriftPosition.hxx>
#include <TUnitsTable.hxx>
#include <THandle.hxx>
#include <TReconHit.hxx>
#include <THit.hxx>
#include <CaptGeomId.hxx>
#include <TRuntimeParameters.hxx>

namespace CP {

  struct hitTimeLowerBound {
        bool operator () (const CP::THandle<CP::THit>& lhs,
                          const double& rhs) {
            return lhs->GetTime() < rhs;
        }
    };

    struct hitTimeUpperBound {
        bool operator () (const double& lhs,
                          const CP::THandle<CP::THit>& rhs) {
            return lhs < rhs->GetTime();
        }
    };

    struct compareHitTimes {
        bool operator () (const CP::THandle<CP::THit>& lhs,
                          const CP::THandle<CP::THit>& rhs) {
            return lhs->GetTime() < rhs->GetTime();
        }
    };

    double hitRMS(CP::THitSelection::iterator begin, 
                  CP::THitSelection::iterator end) {
        double mean = 0;
        double mean2 = 0;
        double w = 0.0;
        while (begin != end) {
            double v = (*begin)->GetCharge();
            mean += v;
            mean2 += v*v;
            w += 1.0;
            ++begin;
        }
        mean /= w;
        mean2 /= w;
        return std::sqrt(mean2 - mean*mean);
    }

    double hitMean(CP::THitSelection::iterator begin, 
                   CP::THitSelection::iterator end) {
        double mean = 0;
        double w = 0.0;
        while (begin != end) {
            double v = (*begin)->GetCharge();
            mean += v;
            w += 1.0;
            ++begin;
        }
        mean /= w;
        return mean;
    }

    double hitTotal(CP::THitSelection::iterator begin, 
                   CP::THitSelection::iterator end) {
        double mean = 0;
        while (begin != end) {
            double v = (*begin)->GetCharge();
            mean += v;
            ++begin;
        }
        return mean;
    }

    /// A base exception for the create track template.
    EXCEPTION(ECreateHit,ECaptRecon);

    /// An exception that there are repeated objects in the track.
    EXCEPTION(EHitRepeatedObject, ECreateHit);

    /// An exception that there are repeated objects in the track.
    EXCEPTION(EHitNonHit, ECreateHit);

  //  template<typename hitIterator>
    bool//CP::THandle<CP::TReconHit> 
    CreateHit(CP::THitSelection& writableHits,
                             const TVector3& hitPosition,
                             double t0,
                             const CP::THandle<CP::THit>& hit1,
                             const CP::THandle<CP::THit>& hit2,
                             const CP::THandle<CP::THit>& hit3);

  void OverlapedXY(const CP::THandle<CP::THit>& hit1,
                               const CP::THandle<CP::THit>& hit2,
                               const CP::THandle<CP::THit>& hit3,
		       TVector3& position);

  const TVector3 PositionXY(const CP::THandle<CP::THit>& hit1,
			    const CP::THandle<CP::THit>& hit2);
};


const double OverlapTime(double r1, double r2, double step){
#ifdef QUADRATURE_OVERLAP
    return std::sqrt(r1*r1 + r2*r2 + step*step);
#else
    return std::max(r1,std::max(r2,step));
#endif
}


const TVector3 CP::PositionXY(const CP::THandle<CP::THit>& hit1,
                                    const CP::THandle<CP::THit>& hit2) {
    double x1 = hit1->GetPosition().X();
    double y1 = hit1->GetPosition().Y();
    double dx1 = hit1->GetYAxis().X();
    double dy1 = hit1->GetYAxis().Y();
    
    double x2 = hit2->GetPosition().X();
    double y2 = hit2->GetPosition().Y();
    double dx2 = hit2->GetYAxis().X();
    double dy2 = hit2->GetYAxis().Y();

   

    // Solve 
    //      x1 + s1*dx1 = x2 + s2*dx2;
    //      y1 + s1*dy1 = y2 + s2*dy2; 
    // for s1, s2
    
    // Only the first shift is used, but here is the full solution (for s1 and
    // s2):
    double s1 = -(dx2*(y1-y2)+dy2*x2-dy2*x1)/(dx2*dy1-dx1*dy2);
    // double s2 = -(dx1*(y1-y2)+dy1*x2-dy1*x1)/(dx2*dy1-dx1*dy2);

    //std::cout<<"x1="<<x1<<"; y1="<<y1<<"; dx1="<<dx1<<"; dy1="<<dy1<<"; x2="<<x2<<"; y2="<<y2<<"; dx2="<<dx1<<"; dy2="<<dy2<<"; s1="<<s1<<std::endl;
    
    return TVector3( (x1+s1*dx1), (y1+s1*dy1), 0.0);
}


const bool OverlappingHits(CP::THandle<CP::THit> h1,
                                     CP::THandle<CP::THit> h2){
#define START_STOP_OVERLAP
#ifdef START_STOP_OVERLAP
    CP::TDriftPosition drift;
    double b1 = drift.GetTime(*h1) + h1->GetTimeStart() - h1->GetTime();
    double b2 = drift.GetTime(*h2) + h2->GetTimeStart() - h2->GetTime();
    double e1 = drift.GetTime(*h1) + h1->GetTimeStop() - h1->GetTime();
    double e2 = drift.GetTime(*h2) + h2->GetTimeStop() - h2->GetTime();
    if (e1 < b2) return false;
    if (e2 < b1) return false;
    return true;
#else
    CP::TDriftPosition drift;
    double t1 = drift.GetTime(*h1);
    double t2 = drift.GetTime(*h2);
    double r1 = fXSeparation*h1->GetTimeRMS();
    double r2 = fUSeparation*h2->GetTimeRMS();
    return (std::abs(t1-t2) < OverlapTime(r1,r2,fMinSeparation));
#endif
}

void CP::OverlapedXY(const CP::THandle<CP::THit>& hit1,
                               const CP::THandle<CP::THit>& hit2,
                               const CP::THandle<CP::THit>& hit3,
                               TVector3& position) {

  
  
    if (CP::GeomId::Captain::GetWirePlane(hit1->GetGeomId())
        == CP::GeomId::Captain::GetWirePlane(hit2->GetGeomId())) {
        CaptError("Hit1 and Hit2 are illegal");
      
    }

    if (CP::GeomId::Captain::GetWirePlane(hit1->GetGeomId())
        == CP::GeomId::Captain::GetWirePlane(hit3->GetGeomId())) {
        CaptError("Hit1 and Hit3 are illegal");
       
    }
    
    if (CP::GeomId::Captain::GetWirePlane(hit2->GetGeomId())
        == CP::GeomId::Captain::GetWirePlane(hit3->GetGeomId())) {
        CaptError("Hit2 and Hit3 are illegal");
    
	}

    TVector3 p1(PositionXY(hit1,hit2));
    TVector3 p2(PositionXY(hit1,hit3));
    double dist = (p2-p1).Mag();

    //if (dist > 2*unit::mm) return false;
    TVector3 p3(PositionXY(hit2,hit3));

    // Find the position for the 3D hit.  Take the average
    // position of the crossing points as the hit position.  It's
    // at Z=0 with a time offset relative to that position.  This
    // should probably be charge weighted, but this is a
    // conservative "average" point.  It's possible that it should
    // be charged weighted, but I suspect not.  That needs to be
    // answered based on fit residuals and pulls.
    position = p1 + p2 + p3;
    position *= 1.0/3.0;
    //  return true;
}



//////////////////////////////////////////////////////////////////
// IMPLEMENTATION
//////////////////////////////////////////////////////////////////

//template<typename hitIterator>
bool//CP::THandle<CP::TReconHit> 
CP::CreateHit(CP::THitSelection& writableHits,
                             const TVector3& hitPosition,
                             double t0,
                             const CP::THandle<CP::THit>& hit1,
                             const CP::THandle<CP::THit>& hit2,
                             const CP::THandle<CP::THit>& hit3) {
  #ifdef DEBUG
  if(hit1){
     std::cout<<"hit1"<<"; X="<<hit1->GetPosition().X()<<"; Y="<<hit1->GetPosition().Y()<<"; Z="<<hit1->GetPosition().Z()<<"; StartTime="<<hit1->GetTimeLowerBound()<<"; StopTime="<<hit1->GetTimeUpperBound()<<std::endl;
     std::cout<<"Wire1#"<<CP::GeomId::Captain::GetWireNumber(hit1->GetGeomId())<<std::endl;
  }
  if(hit2){
     std::cout<<"hit2"<<"; X="<<hit2->GetPosition().X()<<"; Y="<<hit2->GetPosition().Y()<<"; Z="<<hit2->GetPosition().Z()<<"; StartTime="<<hit2->GetTimeLowerBound()<<"; StopTime="<<hit2->GetTimeUpperBound()<<std::endl;
     std::cout<<"Wire2#"<<CP::GeomId::Captain::GetWireNumber(hit2->GetGeomId())<<std::endl;
  }
  if(hit3){
    std::cout<<"hit3"<<"; X="<<hit3->GetPosition().X()<<"; Y="<<hit3->GetPosition().Y()<<"; Z="<<hit3->GetPosition().Z()<<"; StartTime="<<hit3->GetTimeLowerBound()<<"; StopTime="<<hit3->GetTimeUpperBound()<<std::endl;
    std::cout<<"Wire3#"<<CP::GeomId::Captain::GetWireNumber(hit3->GetGeomId())<<std::endl;
  }
  #endif
  /* double fMaxDrift
        = CP::TRuntimeParameters::Get().GetParameterD(
            "captRecon.cluster3d.maxDrift");

    double fXSeparation = CP::TRuntimeParameters::Get().GetParameterD(
            "captRecon.cluster3d.xSeparation");
    double fVSeparation = CP::TRuntimeParameters::Get().GetParameterD(
            "captRecon.cluster3d.vSeparation");
    double fUSeparation = CP::TRuntimeParameters::Get().GetParameterD(
    "captRecon.cluster3d.uSeparation");*/
    double fMinimumOverlap = CP::TRuntimeParameters::Get().GetParameterD(
            "captRecon.cluster3d.minimumOverlap");
    double fMaximumSpread = CP::TRuntimeParameters::Get().GetParameterD(
            "captRecon.cluster3d.maximumSpread");
    
    // double fEnergyPerCharge = CP::TRuntimeParameters::Get().GetParameterD(
    //  "captRecon.energyPerCharge");

    // The time step per digitizer sample.
    double fDigitStep = 500*unit::ns;
    
    // This is the determined by the minimum tick of the digitizer.  
    // double fMinSeparation = fDigitStep;

  

        CP::TDriftPosition drift;
 double startTime = 9E+30;
    double stopTime = -9E+30;
    double expectedCharge = 9E+30;

    bool overlaped=true;
    
    if (hit1 && hit2) {
        if (CP::GeomId::Captain::GetWirePlane(hit1->GetGeomId())
            == CP::GeomId::Captain::GetWirePlane(hit2->GetGeomId())) {
            return false;
        }
        if (hit1->GetTimeLowerBound() > hit2->GetTimeUpperBound()) {
            overlaped= false;
        }
        if (hit1->GetTimeUpperBound() < hit2->GetTimeLowerBound()) {
            overlaped= false;
        }
    }

    if (hit1 && hit3) {
        if (CP::GeomId::Captain::GetWirePlane(hit1->GetGeomId())
            == CP::GeomId::Captain::GetWirePlane(hit3->GetGeomId())) {
            return false;
        }
	 if (hit1->GetTimeLowerBound() > hit3->GetTimeUpperBound()) {
            overlaped= false;
        }
        if (hit1->GetTimeUpperBound() < hit3->GetTimeLowerBound()) {
            overlaped= false;
        }

    }

    if (hit2 && hit3) {
        if (CP::GeomId::Captain::GetWirePlane(hit2->GetGeomId())
            == CP::GeomId::Captain::GetWirePlane(hit3->GetGeomId())) {
            return false;
        }
        if (hit2->GetTimeLowerBound() > hit3->GetTimeUpperBound()) {
            overlaped= false;
        }
        if (hit2->GetTimeUpperBound() < hit3->GetTimeLowerBound()) {
            overlaped= false;
        }
        }
	
if(overlaped){
  //  std::cout<<"PassedTImeCheck"<<std::endl;
    // Find the full range of time covered by this hit.
    if (hit1) {
        startTime = std::min(startTime, hit1->GetTimeLowerBound());
        stopTime = std::max(stopTime, hit1->GetTimeUpperBound());
        expectedCharge = std::min(expectedCharge, hit1->GetCharge());

    }

    if (hit2) {
      startTime = std::min(startTime, hit2->GetTimeLowerBound());
        stopTime = std::max(stopTime, hit2->GetTimeUpperBound());
        expectedCharge = std::min(expectedCharge, hit2->GetCharge());

    }

    if (hit3) {
       startTime = std::min(startTime, hit3->GetTimeLowerBound());
        stopTime = std::max(stopTime, hit3->GetTimeUpperBound());
        expectedCharge = std::min(expectedCharge, hit3->GetCharge());

    }

    // Extend the time range to cover for the drift between the wires.
    startTime -= 15*unit::microsecond;
    stopTime += 15*unit::microsecond;
    


    // Build an "overlap charge distribution"
    int bins = (stopTime-startTime)/fDigitStep + 1;
    // std::cout<<"bins="<<bins<<std::endl;
    std::vector<double> overlap(bins);
  
    if (hit1) {
        double hitStart = drift.GetTime(*hit1);
        hitStart += hit1->GetTimeLowerBound()-hit1->GetTime();
        if (hit1->GetTimeSamples() > 1) {
            int ibin = (hitStart - startTime)/fDigitStep;
          int low =(hit1->GetTimeLowerBound()-hit1->GetTimeStart())/fDigitStep;
	    int up = (hit1->GetTimeUpperBound()-hit1->GetTimeStart())/fDigitStep;
	    // std::cout<<"low1="<<low<<"; up="<<up<<std::endl;
	    int bcount=0;
            for (int i = low; i<up; ++i) {
                // The first hit must exist so just assign it.
	     
                overlap[bcount+ibin] = std::max(0.0,hit1->GetTimeSample(i));
		bcount++;

            }
        }
    }

    if (hit2) {
        double hitStart = drift.GetTime(*hit2);
        hitStart += hit2->GetTimeLowerBound()-hit2->GetTime();
        if (hit2->GetTimeSamples() > 1) {
            int ibin = (hitStart - startTime)/fDigitStep;
	    
            for (int i=0; i<ibin; ++i) overlap[i] = 0.0;

	     int low =(hit2->GetTimeLowerBound()-hit2->GetTimeStart())/fDigitStep;
	    int up = (hit2->GetTimeUpperBound()-hit2->GetTimeStart())/fDigitStep;
	    //std::cout<<"low2="<<low<<"; up="<<up<<std::endl;
	    int bcount=0;
            for (int i = low; i<up; ++i) {
                overlap[bcount+ibin] = std::max(0.0,
                                           std::min(overlap[bcount+ibin],
                                                    hit2->GetTimeSample(i)));
		bcount++;

            }
            for (std::size_t i=ibin+hit2->GetTimeSamples();
                 i<overlap.size(); ++i) {
                overlap[i] = 0.0;
            }
        }
    }

    if (hit3) {
        double hitStart = drift.GetTime(*hit3);
        hitStart += hit3->GetTimeLowerBound()-hit3->GetTime();
        if (hit3->GetTimeSamples() > 1) {
            int ibin = (hitStart - startTime)/fDigitStep;
            for (int i=0; i<ibin; ++i) overlap[i] = 0.0;
	     int low =(hit3->GetTimeLowerBound()-hit3->GetTimeStart())/fDigitStep;
	    int up = (hit3->GetTimeUpperBound()-hit3->GetTimeStart())/fDigitStep;
	    //std::cout<<"low3="<<low<<"; up="<<up<<std::endl;
	    int bcount=0;
            for (int i = low; i<up; ++i) {
                overlap[bcount+ibin] = std::max(0.0,
                                           std::min(overlap[bcount+ibin],
                                                    hit3->GetTimeSample(i)));
		bcount++;

            }
            for (std::size_t i=ibin+hit3->GetTimeSamples();
                 i<overlap.size(); ++i) {
                overlap[i] = 0.0;
            }
        }
    }

    // Find the time of the hit based on the overlaps between the wire hits.
    double hitTime = 0.0;
    double hitTimeRMS = 0.0;
    double hitCharge = 0.0;
    double hitStart = overlap.size();
    double hitStop = 0;
    for (std::size_t i = 0; i<overlap.size(); ++i) {
        double t = startTime + fDigitStep*i;
        hitTime += t*overlap[i];
        hitTimeRMS += t*t*overlap[i];
        hitCharge += overlap[i];
        if (i < hitStart && overlap[i] > 1) hitStart = i;
        if (i > hitStop && overlap[i] > 1) hitStop = i;
    }
    
    // Make sure the hits overlapped by a reasonable amount.
    /*  if (hitCharge < fMinimumOverlap*expectedCharge) {
        CaptNamedInfo("Cluster","Insufficient overlap: "
                     << hitCharge << "/" << expectedCharge
                     << " from "
                     << unit::AsString(startTime+hitStart*fDigitStep, "time")
                     << " to " 
                     << unit::AsString(startTime+hitStop*fDigitStep, "time"));
        CP::TCaptLog::IncreaseIndentation();
        if (hit1) {
            CaptNamedVerbose(
                "Cluster",
                "Hit1: ("
                << CP::GeomId::Captain::GetWirePlane(hit1->GetGeomId())
                << "-" << CP::GeomId::Captain::GetWireNumber(hit1->GetGeomId())
                << ")"
                << " " << unit::AsString(hit1->GetCharge(), "pe")
                << " from " 
                << unit::AsString(hit1->GetTimeStart(), "time")
                << " to " 
                << unit::AsString(hit1->GetTimeStop(), "time"));
        }
        if (hit2) {
            CaptNamedVerbose(
                "Cluster",
                "Hit2: ("
                << CP::GeomId::Captain::GetWirePlane(hit2->GetGeomId())
                << "-" << CP::GeomId::Captain::GetWireNumber(hit2->GetGeomId())
                << ")"
                << " " << unit::AsString(hit2->GetCharge(), "pe")
                << " from " 
                << unit::AsString(hit2->GetTimeStart(), "time")
                << " to " 
                << unit::AsString(hit2->GetTimeStop(), "time"));
        }
        if (hit3) {
            CaptNamedVerbose(
                "Cluster",
                "Hit3: ("
                << CP::GeomId::Captain::GetWirePlane(hit3->GetGeomId())
                << "-" << CP::GeomId::Captain::GetWireNumber(hit3->GetGeomId())
                << ")"
                << " " << unit::AsString(hit3->GetCharge(), "pe")
                << " from " 
                << unit::AsString(hit3->GetTimeStart(), "time")
                << " to " 
                << unit::AsString(hit3->GetTimeStop(), "time"));
        }
        CP::TCaptLog::DecreaseIndentation();
        return false;
    }*/

    //  std::cout<<"passedOverlapCheck"<<std::endl;
    // Figure out where the charge overlap should be split.  The initial guess
    // is the beginning and the end (i.e. no spliting at all).
    std::vector<int> splits;
    splits.push_back(hitStart);
    splits.push_back(hitStop+1);

    // Update the original guess at where to divide the current overlapping
    // charge distribution by finding the internal partition points.
    double hitSpread = fDigitStep*(hitStop-hitStart);
    int nSplits = (int) hitSpread/fMaximumSpread + 1;
    if (nSplits > 1) {
        // Determine the number of samples to put into each "split".  This is
        // a double so that the binning effect is spread over all of the
        // splits (an not just in the last one).
        double split = hitSpread/nSplits/fDigitStep;
        for (double s = hitStart + split; s<hitStop-0.5*split; s += split) {
            splits.push_back((int) (s+0.5) );
        }
    }

    // Make sure that the split points are in order.
    std::sort(splits.begin(), splits.end());

    std::vector<int>::iterator begin = splits.begin();
    std::vector<int>::iterator end = begin+1;
    while (end != splits.end()) {
        // Find the time of the hit based on the overlaps between the wire
        // hits.
        double splitTime = 0.0;
        double splitTimeRMS = 0.0;
        double splitCharge = 0.0;
        for (int i = *begin; i< *end; ++i) {
            double t = startTime + fDigitStep*i;
            splitTime += t*overlap[i];
            splitTimeRMS += t*t*overlap[i];
            splitCharge += overlap[i];
	  
        }
	//	std::cout<<"here"<<std::endl;
	//	std::cout<<"SplitTimeProblem?="<<splitTime<<"; SplitCharge?="<<splitCharge<<std::endl;
        splitTime /= splitCharge;
        splitTimeRMS /= splitCharge;
        splitTimeRMS = splitTimeRMS - splitTime*splitTime;
        if (splitTimeRMS > 0) splitTimeRMS = std::sqrt(splitTimeRMS);

        CP::THandle<CP::TWritableReconHit> hit(new CP::TWritableReconHit(hit1,hit2,hit3));
        hit->SetPosition(hitPosition);
        hit->SetTime(splitTime);
	//	std::cout<<"TimeProblem?="<<hit->GetTime()<<std::endl;
        // Correct for the time zero.
	 hit->SetPosition(drift.GetPosition(*hit,t0).Vect());

        hit->SetTime(t0);
        hit->SetTimeUncertainty(splitTimeRMS/sqrt(3.0));
        hit->SetTimeRMS(splitTimeRMS);
        hit->SetCharge(splitCharge);
	//std::cout<<"here1"<<std::endl;
        // Estimate the charge uncertainty.
        double splitChargeUnc = 0.0;
        if (hit1) {
            double w = hit1->GetChargeUncertainty();
            splitChargeUnc += 1.0/(w*w);
        }
        
        if (hit2) {
            double w = hit2->GetChargeUncertainty();
            splitChargeUnc += 1.0/(w*w);
        }
        
        if (hit3) {
            double w = hit3->GetChargeUncertainty();
            splitChargeUnc += 1.0/(w*w);
        }
        
        if (splitChargeUnc > 0.0) {
            hit->SetChargeUncertainty(std::sqrt(1.0/splitChargeUnc));
        }
	//    std::cout<<"here2"<<std::endl;
        // Find the xyRMS and xyUncertainty.  This is not being done
        // correctly, but this should be an acceptable approximation.  The
        // approximation is that the X rms of the three 2D hits is
        // perpendicular to the wire is an estimate of the hit size.  The Z
        // RMS is calculated based on the time RMS of the wire hits.  The XY
        // rms then has the overlap spacing added in to account for lack of
        // perfect overlaps (this is the constant 2mm squared).
        double xyRMS = 0.0;
        if (hit1) xyRMS += hit1->GetRMS().X()*hit1->GetRMS().X();
        if (hit2) xyRMS += hit2->GetRMS().X()*hit2->GetRMS().X();
        if (hit3) xyRMS += hit3->GetRMS().X()*hit3->GetRMS().X();
        xyRMS /= 3.0;
        xyRMS += xyRMS + 4*unit::mm*unit::mm;
        xyRMS = std::sqrt(xyRMS);
        hit->SetRMS(
            TVector3(xyRMS,xyRMS,drift.GetAverageVelocity()*hit->GetTimeRMS()));
	//  std::cout<<"here3"<<std::endl;
        // For the XY uncertainty, assume a uniform position distribution.
        // For the Z uncertainty, just use the time uncertainty.
        hit->SetUncertainty(
            TVector3(2.0*xyRMS/std::sqrt(12.),2.0*xyRMS/std::sqrt(12.),
                     drift.GetAverageVelocity()*hit->GetTimeUncertainty()));

	 CP::THandle<CP::TReconHit> newHit(new CP::TReconHit(*hit));
	 // std::cout<<"Z="<<newHit->GetPosition().Z()<<std::endl;
	  double ht=(newHit->GetConstituent()->GetTime()+1.6*unit::ms)/(500*unit::ns);
	  // std::cout<<"YesOverlaped"<<"ActualTime="<<ht<<"; X="<<newHit->GetPosition().X()<<"; Y="<<newHit->GetPosition().Y()<<"; Z="<<newHit->GetPosition().Z()<<std::endl;
	 if(!std::isnan(newHit->GetPosition().Z())){
	   writableHits.push_back(newHit);}
	//	std::cout<<"here4"<<std::endl;
        begin = end; ++end;
    }
 }else{
  // std::cout<<"NotOverlaped"<<std::endl;
  double newTime=0;
  double newCharge=0;
  int hitCount=0;
  if(hit1){
    newTime+=hit1->GetTime();
    newCharge+=hit1->GetCharge();
    ++hitCount;
  }
   if(hit2){
    newTime+=hit2->GetTime();
    newCharge+=hit2->GetCharge();
    ++hitCount;
  }
    if(hit3){
    newTime+=hit3->GetTime();
    newCharge+=hit3->GetCharge();
    ++hitCount;
  }
    
    //newTime = hit1->GetTime();
     newTime/=hitCount*1.0;
    newCharge/=hitCount*1.0;

   
  
        CP::THandle<CP::TWritableReconHit> hit(new CP::TWritableReconHit(hit1,hit2,hit3));
        hit->SetPosition(hitPosition);
        hit->SetTime(newTime);

        // Correct for the time zero.
	 hit->SetPosition(drift.GetPosition(*hit,t0).Vect());

        hit->SetTime(t0);
        hit->SetTimeUncertainty(hit1->GetTimeUncertainty());
        hit->SetTimeRMS(hit1->GetTimeRMS());
        hit->SetCharge(newCharge);

        // Estimate the charge uncertainty.
        double splitChargeUnc = 0.0;
        if (hit1) {
            double w = hit1->GetChargeUncertainty();
            splitChargeUnc += 1.0/(w*w);
        }
        
        if (hit2) {
            double w = hit2->GetChargeUncertainty();
            splitChargeUnc += 1.0/(w*w);
        }
        
        if (hit3) {
            double w = hit3->GetChargeUncertainty();
            splitChargeUnc += 1.0/(w*w);
        }
        
        if (splitChargeUnc > 0.0) {
            hit->SetChargeUncertainty(std::sqrt(1.0/splitChargeUnc));
        }
        double xyRMS = 0.0;
        if (hit1) xyRMS += hit1->GetRMS().X()*hit1->GetRMS().X();
        if (hit2) xyRMS += hit2->GetRMS().X()*hit2->GetRMS().X();
        if (hit3) xyRMS += hit3->GetRMS().X()*hit3->GetRMS().X();
        xyRMS /= 3.0;
        xyRMS += xyRMS + 4*unit::mm*unit::mm;
        xyRMS = std::sqrt(xyRMS);
        hit->SetRMS(
            TVector3(xyRMS,xyRMS,drift.GetAverageVelocity()*hit->GetTimeRMS()));
	//  std::cout<<"here3"<<std::endl;
        // For the XY uncertainty, assume a uniform position distribution.
        // For the Z uncertainty, just use the time uncertainty.
        hit->SetUncertainty(
            TVector3(2.0*xyRMS/std::sqrt(12.),2.0*xyRMS/std::sqrt(12.),
                     drift.GetAverageVelocity()*hit->GetTimeUncertainty()));

	 CP::THandle<CP::TReconHit> newHit(new CP::TReconHit(*hit));
	  double ht=(newHit->GetConstituent()->GetTime()+1.6*unit::ms)/(500*unit::ns);
	  // std::cout<<"NotOverlaped"<<"ActualTime="<<ht<<"; X="<<newHit->GetPosition().X()<<"; Y="<<newHit->GetPosition().Y()<<"; Z="<<newHit->GetPosition().Z()<<std::endl;
	 //	 std::cout<<"Z="<<newHit->GetPosition().Z()<<std::endl;
	 if(!std::isnan(newHit->GetPosition().Z())){
	
	      writableHits.push_back(newHit);
	 }

 }
// std::cout<<"Hit"<<std::endl;
    return true;
 
}
#endif
