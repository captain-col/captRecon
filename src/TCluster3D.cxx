#include "TCluster3D.hxx"
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

TVector3 CP::TCluster3D::PositionXY(const CP::THandle<CP::THit>& hit1,
                                    const CP::THandle<CP::THit>& hit2) const {
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
    
    return TVector3( (x1+s1*dx1), (y1+s1*dy1), 0.0);
}

double CP::TCluster3D::ChargeOverlap(const CP::THandle<CP::THit>& hit1,
                                     const CP::THandle<CP::THit>& hit2) const {
    double startTime = 9E+30;
    double stopTime = -9E+30;
    double expectedCharge = 9E+30;
    
    // Find the full range of time covered by this hit.
    if (hit1) {
        startTime = std::min(startTime, hit1->GetTimeStart());
        stopTime = std::max(stopTime, hit1->GetTimeStop());
        expectedCharge = std::min(expectedCharge, hit1->GetCharge());
    }

    if (hit2) {
        startTime = std::min(startTime, hit2->GetTimeStart());
        stopTime = std::max(stopTime, hit2->GetTimeStop());
        expectedCharge = std::min(expectedCharge, hit2->GetCharge());
    }

    // Extend the time range to cover for the drift between the wires.
    startTime -= 15*unit::microsecond;
    stopTime += 15*unit::microsecond;

    // Build an "overlap charge distribution"
    int bins = (stopTime-startTime)/fDigitStep + 1;
    std::vector<double> overlap(bins);

    CP::TDriftPosition drift;
    if (hit1) {
        double hitStart = drift.GetTime(*hit1);
        hitStart += hit1->GetTimeStart()-hit1->GetTime();
        if (hit1->GetTimeSamples() > 1) {
            int ibin = (hitStart - startTime)/fDigitStep;
            for (int i = 0; i<hit1->GetTimeSamples(); ++i) {
                // The first hit must exist so just assign it.
                overlap[i+ibin] = std::max(0.0,hit1->GetTimeSample(i));
            }
        }
    }

    if (hit2) {
        double hitStart = drift.GetTime(*hit2);
        hitStart += hit2->GetTimeStart()-hit2->GetTime();
        if (hit2->GetTimeSamples() > 1) {
            int ibin = (hitStart - startTime)/fDigitStep;
            for (int i=0; i<ibin; ++i) overlap[i] = 0.0;
            for (int i = 0; i<hit2->GetTimeSamples(); ++i) {
                overlap[i+ibin] = std::max(0.0,
                                           std::min(overlap[i+ibin],
                                                    hit2->GetTimeSample(i)));
            }
            for (std::size_t i=ibin+hit2->GetTimeSamples();
                 i<overlap.size(); ++i) {
                overlap[i] = 0.0;
            }
        }
    }

    double charge = 0.0;
    for (std::size_t i = 0; i<overlap.size(); ++i) {
        charge += overlap[i];
    }

    return charge;
}

bool CP::TCluster3D::OverlapXY(const CP::THandle<CP::THit>& hit1,
                               const CP::THandle<CP::THit>& hit2,
                               const CP::THandle<CP::THit>& hit3,
                               TVector3& position) const {
    
    if (CP::GeomId::Captain::GetWirePlane(hit1->GetGeomId())
        == CP::GeomId::Captain::GetWirePlane(hit2->GetGeomId())) {
        CaptError("Hit1 and Hit2 are illegal");
        return false;
    }

    if (CP::GeomId::Captain::GetWirePlane(hit1->GetGeomId())
        == CP::GeomId::Captain::GetWirePlane(hit3->GetGeomId())) {
        CaptError("Hit1 and Hit3 are illegal");
        return false;
    }
    
    if (CP::GeomId::Captain::GetWirePlane(hit2->GetGeomId())
        == CP::GeomId::Captain::GetWirePlane(hit3->GetGeomId())) {
        CaptError("Hit2 and Hit3 are illegal");
        return false;
    }

    TVector3 p1(PositionXY(hit1,hit2));
    TVector3 p2(PositionXY(hit1,hit3));
    double dist = (p2-p1).Mag();

    if (dist > 2*unit::mm) return false;
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
    return true;
}

double CP::TCluster3D::TimeZero(const CP::THitSelection& pmts,
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

CP::TCluster3D::TCluster3D()
    : TAlgorithm("TCluster3D", "Cluster Wire Hits") {
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

CP::TCluster3D::~TCluster3D() { }

double CP::TCluster3D::OverlapTime(double r1, double r2, double step) const {
#ifdef QUADRATURE_OVERLAP
    return std::sqrt(r1*r1 + r2*r2 + step*step);
#else
    return std::max(r1,std::max(r2,step));
#endif
}

bool CP::TCluster3D::OverlappingHits(CP::THandle<CP::THit> h1,
                                     CP::THandle<CP::THit> h2) const {
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

bool CP::TCluster3D::MakeHit(CP::THitSelection& writableHits,
                             const TVector3& hitPosition,
                             double t0,
                             const CP::THandle<CP::THit>& hit1,
                             const CP::THandle<CP::THit>& hit2,
                             const CP::THandle<CP::THit>& hit3) const {

    double startTime = 9E+30;
    double stopTime = -9E+30;
    double expectedCharge = 9E+30;
    
    if (hit1 && hit2) {
        if (CP::GeomId::Captain::GetWirePlane(hit1->GetGeomId())
            == CP::GeomId::Captain::GetWirePlane(hit2->GetGeomId())) {
            return false;
        }
        if (hit1->GetTimeStart() > hit2->GetTimeStop()) {
            return false;
        }
        if (hit1->GetTimeStop() < hit2->GetTimeStart()) {
            return false;
        }
    }

    if (hit1 && hit3) {
        if (CP::GeomId::Captain::GetWirePlane(hit1->GetGeomId())
            == CP::GeomId::Captain::GetWirePlane(hit3->GetGeomId())) {
            return false;
        }
        if (hit1->GetTimeStart() > hit3->GetTimeStop()) {
            return false;
        }
        if (hit1->GetTimeStop() < hit3->GetTimeStart()) {
            return false;
        }
    }

    if (hit2 && hit3) {
        if (CP::GeomId::Captain::GetWirePlane(hit2->GetGeomId())
            == CP::GeomId::Captain::GetWirePlane(hit3->GetGeomId())) {
            return false;
        }
        if (hit2->GetTimeStart() > hit3->GetTimeStop()) {
            return false;
        }
        if (hit2->GetTimeStop() < hit3->GetTimeStart()) {
            return false;
        }
    }

    // Find the full range of time covered by this hit.
    if (hit1) {
        startTime = std::min(startTime, hit1->GetTimeStart());
        stopTime = std::max(stopTime, hit1->GetTimeStop());
        expectedCharge = std::min(expectedCharge, hit1->GetCharge());
    }

    if (hit2) {
        startTime = std::min(startTime, hit2->GetTimeStart());
        stopTime = std::max(stopTime, hit2->GetTimeStop());
        expectedCharge = std::min(expectedCharge, hit2->GetCharge());
    }

    if (hit3) {
        startTime = std::min(startTime, hit3->GetTimeStart());
        stopTime = std::max(stopTime, hit3->GetTimeStop());
        expectedCharge = std::min(expectedCharge, hit3->GetCharge());
    }

    // Extend the time range to cover for the drift between the wires.
    startTime -= 15*unit::microsecond;
    stopTime += 15*unit::microsecond;
    
    CP::TDriftPosition drift;

    // Build an "overlap charge distribution"
    int bins = (stopTime-startTime)/fDigitStep + 1;
    std::vector<double> overlap(bins);

    if (hit1) {
        double hitStart = drift.GetTime(*hit1);
        hitStart += hit1->GetTimeStart()-hit1->GetTime();
        if (hit1->GetTimeSamples() > 1) {
            int ibin = (hitStart - startTime)/fDigitStep;
            for (int i = 0; i<hit1->GetTimeSamples(); ++i) {
                // The first hit must exist so just assign it.
                overlap[i+ibin] = std::max(0.0,hit1->GetTimeSample(i));
            }
        }
    }

    if (hit2) {
        double hitStart = drift.GetTime(*hit2);
        hitStart += hit2->GetTimeStart()-hit2->GetTime();
        if (hit2->GetTimeSamples() > 1) {
            int ibin = (hitStart - startTime)/fDigitStep;
            for (int i=0; i<ibin; ++i) overlap[i] = 0.0;
            for (int i = 0; i<hit2->GetTimeSamples(); ++i) {
                overlap[i+ibin] = std::max(0.0,
                                           std::min(overlap[i+ibin],
                                                    hit2->GetTimeSample(i)));
            }
            for (std::size_t i=ibin+hit2->GetTimeSamples();
                 i<overlap.size(); ++i) {
                overlap[i] = 0.0;
            }
        }
    }

    if (hit3) {
        double hitStart = drift.GetTime(*hit3);
        hitStart += hit3->GetTimeStart()-hit3->GetTime();
        if (hit3->GetTimeSamples() > 1) {
            int ibin = (hitStart - startTime)/fDigitStep;
            for (int i=0; i<ibin; ++i) overlap[i] = 0.0;
            for (int i = 0; i<hit3->GetTimeSamples(); ++i) {
                overlap[i+ibin] = std::max(0.0,
                                           std::min(overlap[i+ibin],
                                                    hit3->GetTimeSample(i)));
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
    if (hitCharge < fMinimumOverlap*expectedCharge) {
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
            new CP::TWritableReconHit(hit1,hit2,hit3));
        hit->SetPosition(hitPosition);
        hit->SetTime(splitTime);
        // Correct for the time zero.
        hit->SetPosition(drift.GetPosition(*hit,t0).Vect());
        hit->SetTime(t0);
        hit->SetTimeUncertainty(splitTimeRMS/sqrt(3.0));
        hit->SetTimeRMS(splitTimeRMS);
        hit->SetCharge(splitCharge);
        
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

double
CP::TCluster3D::FindOverlap(const CP::THandle<CP::THit>& hit,
                            const CP::THandle<CP::THit>& constituent) const {
    double total = 0.0;
    double frac = 0.0;
    CP::TDriftPosition drift;
    double hitTime = drift.GetTime(*hit);
    double startTime = drift.GetTime(*constituent);
    startTime += constituent->GetTimeStart() - constituent->GetTime();
    for (int i = 0; i<constituent->GetTimeSamples(); ++i) {
        total += constituent->GetTimeSample(i);
        double t = startTime + fDigitStep*i - hitTime;
        t /= hit->GetTimeRMS();
        if (t > 2.0) continue;
        frac += std::exp(-0.5*t*t)*constituent->GetTimeSample(i);
    }

    if (total>0.0) return frac/total;
    return 1E-5;
}
                            
CP::THandle<CP::TAlgorithmResult>
CP::TCluster3D::Process(const CP::TAlgorithmResult& wires,
                        const CP::TAlgorithmResult& pmts,
                        const CP::TAlgorithmResult&) {
    CaptLog("TCluster3D Process " << GetEvent().GetContext());
    CP::THandle<CP::THitSelection> wireHits = wires.GetHits();
    if (!wireHits) {
        CaptError("No input hits");
        return CP::THandle<CP::TAlgorithmResult>();
    }

#define CHECK_FOR_OVERLAPPED_HITS
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
    std::unique_ptr<CP::THitSelection> used(new CP::THitSelection("used"));
    std::unique_ptr<CP::THitSelection> unused(new CP::THitSelection("unused"));
    std::unique_ptr<CP::THitSelection> clustered(new CP::THitSelection(
                                                   "clustered"));

    /// Only save the used and unused hits if there are fewer than this many
    /// 2D hits.
    const std::size_t hitLimit = 3000;

    std::vector<float> allRMS;
    CP::THitSelection xHits;
    CP::THitSelection vHits;
    CP::THitSelection uHits;
    for (CP::THitSelection::iterator h2 = wireHits->begin(); 
         h2 != wireHits->end(); ++h2) {
        int plane = CP::GeomId::Captain::GetWirePlane((*h2)->GetGeomId());
        allRMS.push_back((*h2)->GetTimeRMS());
        if (plane == CP::GeomId::Captain::kXPlane) {
            xHits.push_back(*h2);
        }
        else if (plane == CP::GeomId::Captain::kVPlane) {
            vHits.push_back(*h2);
        }
        else if (plane == CP::GeomId::Captain::kUPlane) {
            uHits.push_back(*h2);
        }
        else {
            CaptError("Invalid wire plane");
        }
        if (wireHits->size() < hitLimit) unused->push_back(*h2);
    }

    // Set the maximum time difference between 2D clusters that might become a
    // 3D hit.  This is set to be large (i.e. clusters that are spacially
    // separated by more than 25 mm.
    double maxDeltaT = 16*unit::microsecond;

    // The time range of the search needs to be limited when there are a lot
    // of hits.  The form below works well for large numbers of hits, but
    // needs to be tuned for smaller numbers.  It's removed from the
    // calculation for now (2014/2/16) since I'm working on small events.  The
    // problem is that the full, unoptimized, calculation is approximately
    // O(nHits^3), so for large events this is a very slow.
    double deltaRMS = 3.0;
    std::sort(allRMS.begin(), allRMS.end());
    maxDeltaT = std::min(maxDeltaT, deltaRMS*allRMS[0.90*allRMS.size()]);
    maxDeltaT = std::max(maxDeltaT,2*unit::microsecond);
    
    std::sort(xHits.begin(), xHits.end(), compareHitTimes());
    std::sort(vHits.begin(), vHits.end(), compareHitTimes());
    std::sort(uHits.begin(), uHits.end(), compareHitTimes());

    CP::TCaptLog::IncreaseIndentation();
    CaptLog("X Hits: " << xHits.size()
            << " V Hits: " << vHits.size()
            << " U Hits: " << uHits.size()
            << "  max(RMS): " << unit::AsString(maxDeltaT,"time"));
    CP::TCaptLog::DecreaseIndentation();

    CP::TDriftPosition drift;
    CP::THitSelection xvMatch;
    CP::THitSelection xuMatch;
    
    CP::THitSelection writableHits;
    for (CP::THitSelection::iterator xh=xHits.begin();
         xh!=xHits.end(); ++xh) {
        int hitsForThisXHit = 0;

        xvMatch.clear();
#ifdef LOOK_AT_ALL_HITS
        CP::THitSelection::iterator xvBegin = vHits.begin();
        CP::THitSelection::iterator xvEnd = vHits.end();
#else
        CP::THitSelection::iterator xvBegin
            = std::lower_bound(vHits.begin(), vHits.end(),
                               (*xh)->GetTime()-maxDeltaT,
                               hitTimeLowerBound());
        CP::THitSelection::iterator xvEnd
            = std::upper_bound(xvBegin, vHits.end(),
                               (*xh)->GetTime()+maxDeltaT,
                               hitTimeUpperBound());
#endif
        for (CP::THitSelection::iterator vh=xvBegin;
             vh!=xvEnd; ++vh) {
            // Check that the X and V wires overlap in time.
            if (!OverlappingHits(*xh, *vh)) continue;
            xvMatch.push_back(*vh);
        }
        
        xuMatch.clear();
#ifdef LOOK_AT_ALL_HITS
        CP::THitSelection::iterator xuBegin = uHits.begin();
        CP::THitSelection::iterator xuEnd = uHits.end();
#else
        CP::THitSelection::iterator xuBegin
            = std::lower_bound(uHits.begin(), uHits.end(),
                               (*xh)->GetTime()-maxDeltaT,
                               hitTimeLowerBound());
        CP::THitSelection::iterator xuEnd
            = std::upper_bound(xuBegin, uHits.end(),
                               (*xh)->GetTime()+maxDeltaT,
                               hitTimeUpperBound());
#endif
        for (CP::THitSelection::iterator uh=xuBegin;
             uh!=xuEnd; ++uh) {
            // Check that the X and U wires overlap in time.
            if (!OverlappingHits(*xh, *uh)) continue;
            xuMatch.push_back(*uh);
            
            for (CP::THitSelection::iterator vh=xvMatch.begin();
                 vh != xvMatch.end(); ++vh) {
                // Check that the U and V wires overlap in time.
                if (!OverlappingHits(*uh, *vh)) continue;
            
                // Find the points at which the wires cross and check that the
                // wires all cross and one "point".  Two millimeters is a
                // "magic number chosen based on the geometry for a 3mm
                // separation between the wires.  It needs to change if the
                // wire spacing changes.
                TVector3 hitPosition;
                if (!OverlapXY(*xh,*vh,*uh,hitPosition)) continue;
            
                // Create new writables hits from the overlapping hits.  If
                // the 2D wire hits overlap for a long time (several
                // microseconds), this create several new 3D hits.
                if (!MakeHit(writableHits,hitPosition,t0,*xh,*vh,*uh)) continue;
                ++hitsForThisXHit;
                
                // These three wire hits make a 3D point.  Get them into the
                // correct hit selections.  To protect against bogus events
                // that have a to much "garbage", this only gets done if there
                // aren't to many input hits (i.e. wireHits).  That prevents
                // this class from getting into a "to much memory" situation.
                if (wireHits->size() < hitLimit) {
                    used->AddHit(*xh);
                    unused->RemoveHit(*xh);
                    
                    used->AddHit(*vh);
                    unused->RemoveHit(*vh);
                    
                    used->AddHit(*uh);
                    unused->RemoveHit(*uh);
                }
            }
        }
        
        /// Check to see if there are any hits currently associated with the X
        /// hit.  If there are not, then that means there wasn't a "UVX"
        /// triplet with matching times.  It's still possible that there may
        /// be a good UX, or VX hit.  This makes up for the inefficiency to
        /// find a U or V (induction) hit when the track crosses the wire at
        /// an oblique angle.
        if (hitsForThisXHit>0) continue;
        
        CP::THandle<CP::THit> bestPartner;
        double partnerCharge = 0.0;
        for (CP::THitSelection::iterator ph = xuMatch.begin();
             ph != xuMatch.end(); ++ph) {
            if (!bestPartner) {
                bestPartner = *ph;
                partnerCharge = ChargeOverlap(*xh,*ph);
                continue;
            }
            double charge = ChargeOverlap(*xh,*ph);
            if (partnerCharge > charge) continue;
            partnerCharge = charge;
            bestPartner = *ph;
        }
        
        for (CP::THitSelection::iterator ph = xvMatch.begin();
             ph != xvMatch.end(); ++ph) {
            if (!bestPartner) {
                bestPartner = *ph;
                partnerCharge = ChargeOverlap(*xh,*ph);
                continue;
            }
            double charge = ChargeOverlap(*xh,*ph);
            if (partnerCharge > charge) continue;
            partnerCharge = charge;
            bestPartner = *ph;
        }
        
        // No partner was found.
        if (!bestPartner) continue;
        
        // See if it's a reasonable hit.
        if (!MakeHit(writableHits,PositionXY(*xh,bestPartner),
                     t0,*xh,bestPartner,CP::THandle<CP::THit>())) continue;

        // Now check if it can be added to the bookkeeping objects.
        if (wireHits->size() > hitLimit) continue;
        used->AddHit(*xh);
        unused->RemoveHit(*xh);
                
        used->AddHit(bestPartner);
        unused->RemoveHit(bestPartner);
    }

    CaptNamedLog("Cluster","Number of 3D Hits: " << writableHits.size());

#ifdef REMOVE_OUTLIERS
    CP::TRemoveOutliers outliers;
    outliers.Apply(writableHits);
#endif

#define SHARE_CLUSTERED_CHARGE
#ifdef SHARE_CLUSTERED_CHARGE
    // Share the charge among the 3D hits so that the total charge in the
    // event is not overcounted.
    CP::TDistributeCharge share;

    // Fill the charge sharing object.
    for (CP::THitSelection::iterator h = writableHits.begin();
         h != writableHits.end(); ++h) {
        CP::DistributeCharge::TMeasurementGroup& group = share.AddGroup(*h);
        CP::THandle<CP::TWritableReconHit> groupHit = *h;
        for (int i=0; i<groupHit->GetConstituentCount(); ++i) {
            CP::THandle<CP::THit> hit = groupHit->GetConstituent(i); 
            double physicsWeight = FindOverlap(groupHit,hit);
            group.AddMeasurement(hit, hit->GetCharge(), physicsWeight);
        }
    }

    int iterations = 10 + 100000/writableHits.size();
    iterations = std::min(iterations,5000);
    share.Solve(0.01,iterations);

    // Loops over the measurement groups, and update the charges of the 3D
    // hits.  Since the 3D hit handles reference the hits in the writableHits
    // THitSelection, this also updates the hits that will be copied into the
    // output.
    for (CP::TDistributeCharge::Groups::const_iterator g
             = share.GetGroups().begin();
         g != share.GetGroups().end(); ++g) {
        CP::THandle<CP::TWritableReconHit> groupHit = g->GetObject();
        double totalCharge = g->GetGroupCharge();
        double totalSigma = 0.0;
        for(CP::DistributeCharge::TLinks::const_iterator c
                = g->GetLinks().begin();
            c != g->GetLinks().end(); ++c) {
            CP::THandle<CP::THit> hit = (*c)->GetMeasurement()->GetObject();
            // Notice that the sigma is not reduced by the weight.  This is an
            // attempt to capture some of the extra charge error introduced by
            // the charge sharing, but it's not formally correct.
            double sigma = hit->GetChargeUncertainty();
            totalSigma += 1.0/(sigma*sigma);
        }
        totalSigma = std::sqrt(1.0/totalSigma);
        groupHit->SetCharge(totalCharge);
        groupHit->SetChargeUncertainty(totalSigma);
    }
#endif

    // Copy the selection of writable hits into a selection of recon hits.
    for (CP::THitSelection::iterator h = writableHits.begin();
         h != writableHits.end(); ++h) {
        CP::THandle<CP::TWritableReconHit> hit = *h;
        // Don't include hits that have had all their charge taken away by the
        // charge sharing.  The  10 electron cut corresponds to a hit energy of
        // about 340 eV.
        if (hit->GetCharge() < 10.0) continue;
        clustered->push_back(
            CP::THandle<CP::TReconHit>(new CP::TReconHit(*hit)));
    }

    CaptNamedInfo("Cluster","Mean X Hit Charge " 
            << unit::AsString(hitMean(xHits.begin(), xHits.end()),
                              hitRMS(xHits.begin(), xHits.end()),
                              "pe"));
    CaptNamedInfo("Cluster","Total X Hit Charge " 
            << unit::AsString(hitTotal(xHits.begin(), xHits.end()),
                              "pe"));
    CaptNamedInfo("Cluster","Mean V Hit Charge " 
            << unit::AsString(hitMean(vHits.begin(), vHits.end()),
                              hitRMS(vHits.begin(), vHits.end()),
                              "pe"));
    CaptNamedInfo("Cluster","Total V Hit Charge " 
            << unit::AsString(hitTotal(vHits.begin(), vHits.end()),
                              "pe"));
    CaptNamedInfo("Cluster","Mean U Hit Charge "
            << unit::AsString(hitMean(uHits.begin(), uHits.end()),
                              hitRMS(uHits.begin(), uHits.end()),
                              "pe"));
    CaptNamedInfo("Cluster","Total U Hit Charge " 
            << unit::AsString(hitTotal(uHits.begin(), uHits.end()),
                              "pe"));

    CP::TReconObjectContainer* final = new CP::TReconObjectContainer("final");
    CP::THandle<CP::TReconCluster> usedCluster 
        = CreateCluster("clustered", clustered->begin(), clustered->end());
    final->push_back(usedCluster);
    result->AddResultsContainer(final);
    
    CaptLog("Total hit charge: " 
            << unit::AsString(usedCluster->GetEDeposit(),"pe")
            << " is " 
            << unit::AsString(fEnergyPerCharge*usedCluster->GetEDeposit(),
                              "energy")
            << " from " << clustered->size() << " hits.");
    CaptLog("  Used hits: " << used->size()
            << "  Unused hits: " << unused->size());
    
    if (unused->size() > 0) result->AddHits(unused.release());
    if (used->size() > 0) result->AddHits(used.release());
    result->AddHits(clustered.release());

    return result;
}
