#include "TCluster3D.hxx"
#include "TDriftPosition.hxx"

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
    struct compareHitTime {
        bool operator () (const CP::THandle<CP::THit>& lhs,
                     const CP::THandle<CP::THit>& rhs) {
            return lhs->GetTime() < rhs->GetTime();
        }
    };
};

TVector3 CP::TCluster3D::PositionXY(const CP::THandle<CP::THit>& hit1,
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
    
    return TVector3( (x1+s1*dx1), (y1+s1*dy1), 0.0);
}

double CP::TCluster3D::TimeZero(const CP::THitSelection& pmts,
                                const CP::THitSelection& wires) {

    // The first entry is the weight of the PMT hit.  The second entry is the
    // time that is used to select the time zero.  The hit with the highest
    // weight and earliest time is chosen.
    typedef std::pair<double,double> WeightTime;
    std::vector<WeightTime> weightTimes;
    for (CP::THitSelection::const_iterator p = pmts.begin();
         p != pmts.end(); ++p) {
        WeightTime tw((*p)->GetCharge(), (*p)->GetTime());
        for (CP::THitSelection::const_iterator w = wires.begin();
             w != wires.end(); ++w) {
            if ((*w)->GetTime() < (*p)->GetTime()) continue;
            if ((*w)->GetTime() > (*p)->GetTime() + fMaxDrift) continue;
            tw.first += (*w)->GetCharge();
        }
        weightTimes.push_back(tw);
    }

    std::vector<WeightTime>::iterator m 
        = std::max_element(weightTimes.begin(),weightTimes.end());

    return m->second;
}

CP::TCluster3D::TCluster3D()
    : TAlgorithm("TCluster3D", "Cluster Wire Hits") {
    fMaxDrift
        = CP::TRuntimeParameters::Get().GetParameterD(
            "captRecon.cluster3d.maxDrift");

    // These need to be turned into parameters.
    fXSeparation = 2.0;
    fVSeparation = 2.0;
    fUSeparation = 2.0;
    fMinSeparation = 500*unit::ns;

}

CP::TCluster3D::~TCluster3D() { }

CP::THandle<CP::TAlgorithmResult>
CP::TCluster3D::Process(const CP::TAlgorithmResult& wires,
                        const CP::TAlgorithmResult& pmts,
                        const CP::TAlgorithmResult&) {
    CaptLog("TCluster3D Process " << GetEvent().GetContext());
    CP::THandle<CP::THitSelection> wireHits = wires.GetHitSelection();
    if (!wireHits) {
        CaptError("No input hits");
        return CP::THandle<CP::TAlgorithmResult>();
    }
    CaptLog("Hits in event " << wireHits->size());

    CP::THandle<CP::THitSelection> pmtHits = pmts.GetHitSelection();
    if (!pmtHits) {
        CaptError("No PMT hits provide so time 0 can not be found");
        return CP::THandle<CP::TAlgorithmResult>();
    }

    double t0 = TimeZero(*pmtHits,*wireHits);
    CaptLog("Event time zero is " << t0);
        
    CP::THandle<CP::TAlgorithmResult> result = CreateResult();
    std::auto_ptr<CP::THitSelection> used(new CP::THitSelection("used"));
    std::auto_ptr<CP::THitSelection> unused(new CP::THitSelection("unused"));
    std::auto_ptr<CP::THitSelection> clustered(new CP::THitSelection(
                                                   "clustered"));

    std::auto_ptr<CP::THitSelection> xHits(new CP::THitSelection);
    std::auto_ptr<CP::THitSelection> vHits(new CP::THitSelection);
    std::auto_ptr<CP::THitSelection> uHits(new CP::THitSelection);
    for (CP::THitSelection::iterator h2 = wireHits->begin(); 
         h2 != wireHits->end(); ++h2) {
        int plane = CP::GeomId::Captain::GetWirePlane((*h2)->GetGeomId());
        if (plane == CP::GeomId::Captain::kXPlane) {
            xHits->push_back(*h2);
        }
        else if (plane == CP::GeomId::Captain::kVPlane) {
            vHits->push_back(*h2);
        }
        else if (plane == CP::GeomId::Captain::kUPlane) {
            uHits->push_back(*h2);
        }
        else {
            CaptError("Invalid wire plane");
        }
        unused->push_back(*h2);
    }

    std::sort(xHits->begin(), xHits->end(), compareHitTime());
    std::sort(vHits->begin(), vHits->end(), compareHitTime());
    std::sort(uHits->begin(), uHits->end(), compareHitTime());

    CP::TDriftPosition drift;

    // The time must be drift corrected!
    double deltaT;
    
    for (CP::THitSelection::iterator xh=xHits->begin();
         xh!=xHits->end(); ++xh) {
        double xTime = drift.GetTime(**xh);
        double xRMS = (*xh)->GetTimeRMS();
        for (CP::THitSelection::iterator vh=vHits->begin(); 
             vh!=vHits->end(); ++vh) {
            double vTime = drift.GetTime(**vh);
            deltaT = std::abs(vTime - xTime);
            if (deltaT > 10*unit::microsecond) continue;
            double vRMS = (*vh)->GetTimeRMS();

            for (CP::THitSelection::iterator uh=uHits->begin(); 
                 uh!=uHits->end(); ++uh) {
                double uTime = drift.GetTime(**uh);
                deltaT = std::abs(uTime - xTime);
                if (deltaT > 10*unit::microsecond) continue;
                double uRMS = (*uh)->GetTimeRMS();

                // Check that the X and U wires overlap in time.
                if (deltaT > (fXSeparation*xRMS
                              + fUSeparation*uRMS
                              + fMinSeparation)) continue;
                
                // Check that the X and V wires overlap in time.
                deltaT = std::abs(vTime - xTime);                
                if (deltaT > (fXSeparation*xRMS
                              + fVSeparation*vRMS
                              + fMinSeparation)) continue;

                // Check that the U and V wires overlap in time.
                deltaT = std::abs(vTime - uTime);                
                if (deltaT > (fUSeparation*uRMS
                              + fVSeparation*vRMS
                              + fMinSeparation)) continue;

                // Find the points at which the wires cross and check that the
                // wires all cross and one "point".  Two millimeters is a
                // "magic number chosen based on the geometry for a 3mm
                // separation between the wires.  It needs to change if the
                // wire spacing changes.
                TVector3 p1(PositionXY(*xh,*vh));
                TVector3 p2(PositionXY(*xh,*uh));
                double dist = (p2-p1).Mag();
                if (dist > 2*unit::mm) continue;

                CP::TWritableReconHit hit(*xh,*vh,*uh);
                TVector3 p3(PositionXY(*vh,*uh));

                // These three wire hits make a 3D point.a
                used->AddHit(*xh);
                unused->RemoveHit(*xh);

                used->AddHit(*vh);
                unused->RemoveHit(*vh);

                used->AddHit(*uh);
                unused->RemoveHit(*uh);

                // Set the time.  It's the average of the wire hit times after
                // drifting to Z = 0.0.
                double t1 = drift.GetTime(**xh);
                double t2 = drift.GetTime(**vh);
                double t3 = drift.GetTime(**uh);
                double tAvg = (t1+t2+t3)/3.0;
                hit.SetTime(tAvg); // This will be overridden to be time zero.

                // Find the RMS of the hit time.  This takes into account the
                // RMS of the individual hits as well as the offset of the
                // hits from the average of their times.
                double tRMS = (t1-tAvg)*(t1-tAvg);
                tRMS += (*xh)->GetTimeRMS()*(*xh)->GetTimeRMS();
                tRMS += (t2-tAvg)*(t2-tAvg);
                tRMS += (*vh)->GetTimeRMS()*(*vh)->GetTimeRMS();
                tRMS += (t3-tAvg)*(t3-tAvg);
                tRMS += (*uh)->GetTimeRMS()*(*uh)->GetTimeRMS();
                tRMS /= 3.0;
                tRMS = std::sqrt(tRMS);
                hit.SetTimeRMS(tRMS);
                
                // This isn't exactly right, but it's a pretty good estimate.
                double tUnc = (t1-tAvg)*(t1-tAvg);
                tUnc += (t2-tAvg)*(t2-tAvg);
                tUnc += (t3-tAvg)*(t3-tAvg);
                tUnc /= 3.0;
                tUnc = std::sqrt(tUnc);
                hit.SetTimeUncertainty(tUnc);

                // For now set the charge to X wire charge.  This doesn't work
                // for overlapping hits but gives a reasonable estimate of the
                // energy deposition otherwise.  The U and V wires don't
                // measure the total charge very well.
                hit.SetCharge((*xh)->GetCharge());
                hit.SetChargeUncertainty((*xh)->GetChargeUncertainty());
                
                // Find the position for the 3D hit.  Take the average
                // position of the crossing points as the hit position.  It's
                // at Z=0 with a time offset relative to that position.  This
                // should probably be charge weighted, but this is a
                // conservative "average" point.  It's possible that it should
                // be charged weighted, but I suspect not.  That needs to be
                // answered based on fit residuals and pulls.
                TVector3 pos = p1 + p2 + p3;
                pos *= 1.0/3.0;
                pos.SetZ(0.0);
                hit.SetPosition(pos);

                // Find the xyRMS and xyUncertainty.  This is not being done
                // correctly, but this should be an acceptable approximation
                // for now.  The approximation is based on the idea that the X
                // rms of the three 2D hits is perpendicular to the wire, and
                // takes that as an estimate of the hit size.  The Z RMS is
                // calculated based on the time RMS of the three hits.
                double xyRMS = 0.0;
                xyRMS += (*xh)->GetRMS().X()*(*xh)->GetRMS().X();
                xyRMS += (*vh)->GetRMS().X()*(*vh)->GetRMS().X();
                xyRMS += (*uh)->GetRMS().X()*(*uh)->GetRMS().X();
                xyRMS /= 3.0;
                xyRMS += xyRMS + dist*dist;
                xyRMS = std::sqrt(xyRMS);

                hit.SetRMS(
                    TVector3(xyRMS,xyRMS,
                             drift.GetAverageDriftVelocity()*tRMS));

                // For the XY uncertainty, assume a uniform position
                // distribution.  For the Z uncertainty, just use the time
                // uncertainty.
                hit.SetUncertainty(
                    TVector3(2.0*xyRMS/std::sqrt(12.),2.0*xyRMS/std::sqrt(12.),
                             drift.GetAverageDriftVelocity()*tUnc));

                // Correct for the time zero.
                TLorentzVector driftedPosition = drift.GetPosition(hit,t0);
                hit.SetTime(t0);
                hit.SetPosition(driftedPosition.Vect());

                clustered->push_back(CP::THandle<CP::TReconHit>(
                                         new CP::TReconHit(hit)));
            }
        }
    }
        
    CaptLog("Used Hits: " << used->size());
    CaptLog("Unused Hits: " << unused->size());
    CaptLog("Clustered Hits: " << clustered->size());

    CP::TReconObjectContainer* final = new CP::TReconObjectContainer("final");
    CP::THandle<CP::TReconCluster> usedCluster(new CP::TReconCluster);
    usedCluster->FillFromHits("clustered",*clustered);
    final->push_back(usedCluster);
    result->AddResultsContainer(final);

    result->AddHitSelection(unused.release());
    result->AddHitSelection(used.release());
    result->AddHitSelection(clustered.release());

    return result;
}
