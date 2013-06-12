#include <TCluster3D.hxx>
#include <TDriftPosition.hxx>

#include <THandle.hxx>
#include <TReconHit.hxx>
#include <TReconCluster.hxx>
#include <TCaptLog.hxx>
#include <CaptGeomId.hxx>
#include <HEPUnits.hxx>

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
    
    // The full solution:
    double s1 = -(dx2*(y1-y2)+dy2*x2-dy2*x1)/(dx2*dy1-dx1*dy2);
    // double s2 = -(dx1*(y1-y2)+dy1*x2-dy1*x1)/(dx2*dy1-dx1*dy2); 
    
    return TVector3(x1+s1*dx1,y1+s1*dy1,0.0);
}

CP::TCluster3D::TCluster3D()
    : TAlgorithm("TCluster3D", "Cluster Wire Hits") { }

CP::TCluster3D::~TCluster3D() { }

CP::THandle<CP::TAlgorithmResult>
CP::TCluster3D::Process(const CP::TAlgorithmResult& input) {
    CaptLog("Process " << GetEvent().GetContext());
    CP::THandle<CP::THitSelection> hits2d = input.GetHitSelection();
    if (!hits2d) return CP::THandle<CP::TAlgorithmResult>();

    CP::THandle<CP::TAlgorithmResult> result = CreateResult();
    std::auto_ptr<CP::THitSelection> used(new CP::THitSelection("used"));
    std::auto_ptr<CP::THitSelection> unused(new CP::THitSelection("unused"));
    std::auto_ptr<CP::THitSelection> clustered(new CP::THitSelection(
                                                   "clustered"));

    std::auto_ptr<CP::THitSelection> xHits(new CP::THitSelection);
    std::auto_ptr<CP::THitSelection> vHits(new CP::THitSelection);
    std::auto_ptr<CP::THitSelection> uHits(new CP::THitSelection);
    for (CP::THitSelection::iterator h2 = hits2d->begin(); 
         h2 != hits2d->end(); ++h2) {
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
    int count = 0;
    for (CP::THitSelection::iterator xh=xHits->begin();
         xh!=xHits->end(); ++xh) {
        for (CP::THitSelection::iterator vh=vHits->begin(); 
             vh!=vHits->end(); ++vh) {
            deltaT = std::abs(drift.GetTime(*xh) - drift.GetTime(*vh));
            if (deltaT > 1*unit::microsecond) continue;
            for (CP::THitSelection::iterator uh=uHits->begin(); 
                 uh!=uHits->end(); ++uh) {
                deltaT = std::abs(drift.GetTime(*xh) - drift.GetTime(*uh));
                if (deltaT > 1*unit::microsecond) continue;
                deltaT = std::abs(drift.GetTime(*vh) - drift.GetTime(*uh));
                if (deltaT > 1*unit::microsecond) continue;

                /// Find the points at which the wires cross.
                TVector3 p1(PositionXY(*xh,*vh));
                TVector3 p2(PositionXY(*xh,*uh));
                double dist = (p2-p1).Mag();
                if (dist > 2*unit::mm) continue;

                CP::TWritableReconHit hit(*xh,*vh,*uh);
                TVector3 p3(PositionXY(*vh,*uh));

                // These three wire hits make a 3D point.
                used->AddHit(*xh);
                unused->RemoveHit(*xh);

                used->AddHit(*vh);
                unused->RemoveHit(*vh);

                used->AddHit(*uh);
                unused->RemoveHit(*uh);

                // Take the average position of the crossing points as the hit
                // position.
                TVector3 pos = p1 + p2 + p3;
                pos *= 1/3.0;
                pos.SetZ(0.0);
                hit.SetPosition(pos);

                // Set the RMS and Uncertainty.  This is not being done
                // correctly, and should depend on the time RMS of the input
                // hits.
                hit.SetRMS(TVector3(dist,dist,dist));
                hit.SetUncertainty(hit.GetRMS());

                /// Set the time.
                double t1 = drift.GetTime(*xh);
                double t2 = drift.GetTime(*vh);
                double t3 = drift.GetTime(*uh);

                hit.SetTime((t1+t2+t3)/3.0);
                hit.SetTimeUncertainty(500*unit::ns);
                hit.SetTimeRMS((*xh)->GetTimeRMS());

                hit.SetCharge((*xh)->GetCharge() 
                              + (*vh)->GetCharge()
                              + (*uh)->GetCharge());

                hit.SetChargeUncertainty(
                    std::max((*xh)->GetChargeUncertainty(),
                             std::max((*vh)->GetChargeUncertainty(),
                                      (*uh)->GetChargeUncertainty())));

                // CaptLog("Add " << hit.GetPosition().X()
                //          << " " << hit.GetPosition().Y()
                //          << " " << hit.GetPosition().Z()
                //          << " " << hit.GetTime()
                //          << " " << hit.GetCharge());

                /// Find the position for the 3D hit.  It's at Z=0 with a time
                /// offset relative to that position.
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
