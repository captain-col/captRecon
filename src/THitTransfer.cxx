#include "THitTransfer.hxx"
#include "TDriftPosition.hxx"
#include "CreateCluster.hxx"
#include "TShareCharge.hxx"
#include "TDistributeCharge.hxx"
#include "TRemoveOutliers.hxx"

#include <THandle.hxx>
#include <TReconHit.hxx>
#include <TReconCluster.hxx>
#include <TCaptLog.hxx>
#include <CaptGeomId.hxx>
#include <HEPUnits.hxx>
#include <TUnitsTable.hxx> 
#include <TRuntimeParameters.hxx>

#include <algorithm>
#include <memory>
#include <cmath>

namespace {
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
};


double CP::THitTransfer::TimeZero(const CP::THitSelection& pmts,
                                const CP::THitSelection& wires) {
    ///////////////////////////////////////////////////////////////////////
    // Find the event time zero.
    ///////////////////////////////////////////////////////////////////////

    if (pmts.empty()) return 0.0;
    
    std::vector<double> times;
    for (CP::THitSelection::const_iterator p = pmts.begin();
         p != pmts.end(); ++p) {
        times.push_back((*p)->GetTime());
    }
    std::sort(times.begin(),times.end());

    double t0 = 0.0;
    int maxHits = 0;
    for (std::vector<double>::iterator t = times.begin();
         t != times.end(); ++t) {
        std::vector<double>::iterator h = t;
        while (++h != times.end()) {
            if (*h - *t > 2*unit::microsecond) break;
        }
        int dh = h - t;
        if (dh > maxHits) {
            maxHits = dh;
            t0 = *t;
        }
    }

    return t0;
}

CP::THitTransfer::THitTransfer()
    : TAlgorithm("THitTransfer", "Cluster Wire Hits") {
    fMaxDrift
        = CP::TRuntimeParameters::Get().GetParameterD(
            "captRecon.cluster3d.maxDrift");

    fXSeparation = CP::TRuntimeParameters::Get().GetParameterD(
            "captRecon.cluster3d.xSeparation");
    fVSeparation = CP::TRuntimeParameters::Get().GetParameterD(
            "captRecon.cluster3d.vSeparation");
    fUSeparation = CP::TRuntimeParameters::Get().GetParameterD(
            "captRecon.cluster3d.uSeparation");
    fMinimumOverlap = CP::TRuntimeParameters::Get().GetParameterD(
            "captRecon.cluster3d.minimumOverlap");
    fMaximumSpread = CP::TRuntimeParameters::Get().GetParameterD(
            "captRecon.cluster3d.maximumSpread");   
    fEnergyPerCharge = CP::TRuntimeParameters::Get().GetParameterD(
        "captRecon.energyPerCharge");
    // The time step per digitizer sample.
    fDigitStep = 500*unit::ns;   
    // This is the determined by the minimum tick of the digitizer.  
    fMinSeparation = fDigitStep;

}

CP::THitTransfer::~THitTransfer() { }

bool CP::THitTransfer::MakeHit(CP::THitSelection& writableHits,
                             const TVector3& hitPosition,
                             double t0,
                             const CP::THandle<CP::THit>& hit1) const {

    double startTime = 9E+30;
    double stopTime = -9E+30;
    double expectedCharge = 9E+30;

    
    // std::cout<<"HITTIMEUNS="<<hit1->GetTimeUncertainty()<<std::endl;
    // Find the full range of time covered by this hit.
    if (hit1) {
        startTime = std::min(startTime, hit1->GetTimeLowerBound());
        stopTime = std::max(stopTime, hit1->GetTimeUpperBound());
        expectedCharge = std::min(expectedCharge, hit1->GetCharge());
    }

  
    // Extend the time range to cover for the drift between the wires.
    startTime -= 15*unit::microsecond;
    stopTime += 15*unit::microsecond;
    
    CP::TDriftPosition drift;

    // Build an "overlap charge distribution"
    int bins = (stopTime-startTime)/fDigitStep + 2;

    std::vector<double> overlap(bins);

    if (hit1) {
  
        double hitStart = drift.GetTime(*hit1);

        hitStart += hit1->GetTimeLowerBound()-hit1->GetTime();

        if (hit1->GetTimeSamples() > 1) {

            int ibin = (hitStart - startTime)/fDigitStep;
	    int low =(hit1->GetTimeLowerBound()-hit1->GetTimeStart())/fDigitStep;
	    int up = (hit1->GetTimeUpperBound()-hit1->GetTimeStart())/fDigitStep;

	    int bcount=0;
            for (int i = low; i<up; ++i) {
                // The first hit must exist so just assign it.
	       overlap[bcount+ibin] = std::max(0.0,hit1->GetTimeSample(i));
	       bcount++;
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
        splitTime /= splitCharge;
        splitTimeRMS /= splitCharge;
        splitTimeRMS = splitTimeRMS - splitTime*splitTime;
        if (splitTimeRMS > 0) splitTimeRMS = std::sqrt(splitTimeRMS);
        
        CP::THandle<CP::TWritableReconHit> hit(
            new CP::TWritableReconHit(hit1));
        hit->SetPosition(hitPosition);
        hit->SetTime(splitTime);

        // Correct for the time zero.
        hit->SetPosition(drift.GetPosition(*hit,t0).Vect());
        //hit->SetTime(t0);
        
	       hit->SetTimeUncertainty(splitTimeRMS/sqrt(3.0));
	 hit->SetTimeRMS(splitTimeRMS);
	 hit->SetCharge(splitCharge);
        
	   hit->SetTimeUncertainty(hit1->GetTimeUncertainty());
	   hit->SetTimeRMS(hit1->GetTimeRMS());
	   // hit->SetCharge(0.1);
	   
        
        // Estimate the charge uncertainty.
        double splitChargeUnc = 0.0;
        if (hit1) {
            double w = hit1->GetChargeUncertainty();
            splitChargeUnc += 1.0/(w*w);
        }
        
        if (splitChargeUnc > 0.0) {
            hit->SetChargeUncertainty(std::sqrt(1.0/splitChargeUnc));
        }
        
        // Find the xyRMS and xyUncertainty.  This is not being done
        // correctly, but this should be an acceptable approximation.  The
        // approximation is that the X rms of the three 2D hits is
        // perpendicular to the wire is an estimate of the hit size.  The Z
        // RMS is calculated based on the time RMS of the wire hits.  The XY
        // rms then has the overlap spacing added in to account for lack of
        // perfect overlaps (this is the constant 2mm squared).
        double xyRMS = 0.0;
        if (hit1) xyRMS += hit1->GetRMS().X()*hit1->GetRMS().X();
        xyRMS /= 3.0;
        xyRMS += xyRMS + 4*unit::mm*unit::mm;
        xyRMS = std::sqrt(xyRMS);
        hit->SetRMS(
            TVector3(xyRMS,xyRMS,drift.GetAverageVelocity()*hit->GetTimeRMS()));
        
        // For the XY uncertainty, assume a uniform position distribution.
        // For the Z uncertainty, just use the time uncertainty.
        hit->SetUncertainty(
            TVector3(2.0*xyRMS/std::sqrt(12.),2.0*xyRMS/std::sqrt(12.),
                     drift.GetAverageVelocity()*hit->GetTimeUncertainty()));
        
        writableHits.push_back(hit);
        begin = end; ++end;
    }
    
    return true;
}


                            
CP::THandle<CP::TAlgorithmResult>
CP::THitTransfer::Process(const CP::TAlgorithmResult& wires,
                        const CP::TAlgorithmResult& pmts,
                        const CP::TAlgorithmResult&) {
    CaptLog("THitTransfer Process " << GetEvent().GetContext());
   CP::THandle<CP::THitSelection> wireHits_all = wires.GetHits();
    if (!wireHits_all) {
        CaptError("No input hits");
        return CP::THandle<CP::TAlgorithmResult>();
    }

    CP::THitSelection xHits;
    CP::THitSelection vHits;
    CP::THitSelection uHits;
     CP::THitSelection xHits_all;
    CP::THitSelection vHits_all;
    CP::THitSelection uHits_all;
    for (CP::THitSelection::iterator h = wireHits_all->begin(); 
         h != wireHits_all->end(); ++h) {
        int plane = CP::GeomId::Captain::GetWirePlane((*h)->GetGeomId());
        if (plane == CP::GeomId::Captain::kXPlane) {
	   xHits_all.push_back(*h);
        }
        else if (plane == CP::GeomId::Captain::kVPlane) {
	  vHits_all.push_back(*h);
        }
        else if (plane == CP::GeomId::Captain::kUPlane) {
	  uHits_all.push_back(*h);
        }
    }

         for(CP::THitSelection::iterator i = xHits_all.begin(); i !=xHits_all.end(); i++){
	double SampleValue_X = 0.0;
	for ( int count=0; count < ( (*i)->GetTimeSamples()); count++)
	  { 
	    SampleValue_X = SampleValue_X+((*i)->GetTimeSample(count));
	              }
		double SampleValueAbs_X = 0.0;
	for ( int count=0; count < ( (*i)->GetTimeSamples()); count++)
	  { 
	    SampleValueAbs_X = SampleValueAbs_X+fabs(((*i)->GetTimeSample(count)));
	              }
	if(fabs(SampleValue_X)/SampleValueAbs_X>0.2){
	  xHits.push_back(*i);
	}
           }
    
     for(CP::THitSelection::iterator i = uHits_all.begin(); i !=uHits_all.end(); i++){
	  uHits.push_back(*i);
           }

      for(CP::THitSelection::iterator i = vHits_all.begin(); i !=vHits_all.end(); i++){
        vHits.push_back(*i);  
	 }

      CP::THandle<CP::THitSelection> wireHits(new CP::THitSelection);
      for(CP::THitSelection::iterator i = xHits.begin(); i !=xHits.end(); i++){
	wireHits->push_back(*i);
      }
      for(CP::THitSelection::iterator i = uHits.begin(); i !=uHits.end(); i++){
	wireHits->push_back(*i);
      }
      for(CP::THitSelection::iterator i = vHits.begin(); i !=vHits.end(); i++){
	wireHits->push_back(*i);
      }



    

    //#define CHECK_FOR_OVERLAPPED_HITS
#ifdef  CHECK_FOR_OVERLAPPED_HITS
    for (CP::THitSelection::iterator h = wireHits->begin();
         h != wireHits->end(); ++h) {
        for (CP::THitSelection::iterator i = h+1; i != wireHits->end(); ++i) {
            // The hits aren't the same channel, so they can't overlap.
            if ((*i)->GetChannelId() != (*h)->GetChannelId()) continue;
            // "i" starts after "h" ends, so they can't overlap.
            if ((*i)->GetTimeStart() >= (*h)->GetTimeStop()) continue;
            // "i" ends before "h" begins, so they can't overlap.
            if ((*i)->GetTimeStop() <= (*h)->GetTimeStart()) continue;
            CaptError("Overlapping hit on " << (*h)->GetChannelId());
            CaptError("   hit " << (*h)->GetTimeStart()
                      << " " << (*h)->GetTimeStop());
            CaptError("   hit " << (*i)->GetTimeStart()
                      << " " << (*i)->GetTimeStop());
        }
    }
#endif
   
    double t0 = 0.0;
     CP::THandle<CP::THitSelection>  pmtHits = pmts.GetHits();
     if (pmtHits) t0 = TimeZero(*pmtHits,*wireHits);
      
    CP::THandle<CP::TAlgorithmResult> result = CreateResult();

    typedef std::set< CP::THandle<CP::THit> >  HitSet;
    HitSet usedSet;

    
    CP::THitSelection writableHits;


    // All of the 3 wire 3D hits have been found, now check any unassociated
    // hits to see if there are any 2 wire 3D hits that should be formed.
    for (CP::THitSelection::iterator h = wireHits->begin(); 
         h != wireHits->end(); ++h) {
      TVector3 pos((*h)->GetPosition().X(),(*h)->GetPosition().Y(),0.0);
       int plane = CP::GeomId::Captain::GetWirePlane((*h)->GetGeomId());

       if (plane == CP::GeomId::Captain::kVPlane) {
	 // pos.RotateZ(30);
        }
        else if (plane == CP::GeomId::Captain::kUPlane) {
	  // pos.RotateZ(30);
        }
            if (!MakeHit(writableHits,pos,
                         t0,*h)) continue;
            usedSet.insert(*h);
      

    }

    std::unique_ptr<CP::THitSelection> clustered(new CP::THitSelection(
                                                     "clustered"));
     std::unique_ptr<CP::TReconObjectContainer> final(
        new CP::TReconObjectContainer("final"));
    
    // Copy the selection of writable hits into a selection of recon hits.
    for (CP::THitSelection::iterator h = writableHits.begin();
         h != writableHits.end(); ++h) {
        CP::THandle<CP::TWritableReconHit> hit = *h;
        // Don't include hits that have had all their charge taken away by the
        // charge sharing.  The  10 electron cut corresponds to a hit energy of
        // about 340 eV.
        if (hit->GetCharge() < 10.0) continue;
        CP::THandle<CP::TReconHit> newHit(new CP::TReconHit(*hit));
        clustered->push_back(newHit);
	//	final->push_back(newHit);
       	int plane = CP::GeomId::Captain::GetWirePlane((*newHit).GetGeomId());
	/*	CaptNamedLog("THitTransfering",
		"plane " << plane << " PositionX= " << (*newHit).GetPosition().X() << " PositionY="<<(*newHit).GetPosition().Y()<<" PositionZ="<<(*newHit).GetPosition().Z());*/
      
    }


   
   
    /// Save the used 2D wire hits for output.
    std::unique_ptr<CP::THitSelection> used(new CP::THitSelection("used"));
    for (HitSet::iterator h = usedSet.begin(); h != usedSet.end(); ++h) {
        used->push_back(*h);
    }

    /// Save the unused 2D wire hits for output.
    std::unique_ptr<CP::THitSelection> unused(new CP::THitSelection("unused"));
    for (CP::THitSelection::iterator h = wireHits->begin(); 
         h != wireHits->end(); ++h) {
        if (usedSet.find(*h) != usedSet.end()) continue;
        unused->push_back(*h);
    }

  
    // if (unused->size() > 0) result->AddHits(unused.release());
    // if (used->size() > 0) result->AddHits(used.release());
     if (clustered->size() > 0) result->AddHits(clustered.release());
    result->AddResultsContainer(final.release());

    return result;
}
