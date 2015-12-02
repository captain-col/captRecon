#include <TShareCharge.hxx>

#include <HEPUnits.hxx>
#include <TCaptLog.hxx>
#include <TFADCHit.hxx>
#include <TReconHit.hxx>
#include <CaptGeomId.hxx>

#include <tut.h>

namespace tut {
    struct baseShareCharge {
        baseShareCharge() {
            // Run before each test.
        }
        ~baseShareCharge() {
            // Run after each test.
        }
    };

    // Declare the test
    typedef test_group<baseShareCharge>::object testShareCharge;
    test_group<baseShareCharge> groupShareCharge("TShareCharge");

    // Test the declaration.
    template<> template<> void testShareCharge::test<1> () {
        CP::TShareCharge share;
    }

    // Test the declaration.
    template<> template<> void testShareCharge::test<2> () {
        CP::THitSelection uHits;
        CP::THitSelection vHits;
        CP::THitSelection xHits;

        CP::TWritableFADCHit hit;

        // Make the X hits
        hit.SetGeomId(CP::GeomId::Captain::Wire(
                          CP::GeomId::Captain::kXPlane,1));
        hit.SetCharge(5.0);
        hit.SetChargeUncertainty(1.0);
        hit.SetTime(0.0);
        hit.SetTimeRMS(1.0);
        hit.SetTimeStart(-1.0);
        hit.SetTimeStop(1.0);
        hit.SetTimeUncertainty(1.0);
        xHits.push_back(CP::THandle<CP::TFADCHit>(new CP::TFADCHit(hit)));

        hit.SetGeomId(CP::GeomId::Captain::Wire(
                          CP::GeomId::Captain::kXPlane,2));
        hit.SetCharge(3.0);
        hit.SetChargeUncertainty(1.0);
        hit.SetTime(0.0);
        hit.SetTimeRMS(1.0);
        hit.SetTimeStart(-1.0);
        hit.SetTimeStop(1.0);
        hit.SetTimeUncertainty(1.0);
        xHits.push_back(CP::THandle<CP::TFADCHit>(new CP::TFADCHit(hit)));
        
        // Make the U hits
        hit.SetGeomId(CP::GeomId::Captain::Wire(
                          CP::GeomId::Captain::kUPlane,1));
        hit.SetCharge(5.0);
        hit.SetChargeUncertainty(1.0);
        hit.SetTime(0.0);
        hit.SetTimeRMS(1.0);
        hit.SetTimeStart(-1.0);
        hit.SetTimeStop(1.0);
        hit.SetTimeUncertainty(1.0);
        uHits.push_back(CP::THandle<CP::TFADCHit>(new CP::TFADCHit(hit)));

        hit.SetGeomId(CP::GeomId::Captain::Wire(
                          CP::GeomId::Captain::kUPlane,2));
        hit.SetCharge(3.0);
        hit.SetChargeUncertainty(1.0);
        hit.SetTime(0.0);
        hit.SetTimeRMS(1.0);
        hit.SetTimeStart(-1.0);
        hit.SetTimeStop(1.0);
        hit.SetTimeUncertainty(1.0);
        uHits.push_back(CP::THandle<CP::TFADCHit>(new CP::TFADCHit(hit)));
        
        // Make the V hits
        hit.SetGeomId(CP::GeomId::Captain::Wire(
                          CP::GeomId::Captain::kVPlane,1));
        hit.SetCharge(7.0);
        hit.SetChargeUncertainty(1.0);
        hit.SetTime(0.0);
        hit.SetTimeRMS(1.0);
        hit.SetTimeStart(-1.0);
        hit.SetTimeStop(1.0);
        hit.SetTimeUncertainty(1.0);
        vHits.push_back(CP::THandle<CP::TFADCHit>(new CP::TFADCHit(hit)));

        hit.SetGeomId(CP::GeomId::Captain::Wire(
                          CP::GeomId::Captain::kVPlane,2));
        hit.SetCharge(1.0);
        hit.SetChargeUncertainty(1.0);
        hit.SetTime(0.0);
        hit.SetTimeRMS(1.0);
        hit.SetTimeStart(-1.0);
        hit.SetTimeStop(1.0);
        hit.SetTimeUncertainty(1.0);
        vHits.push_back(CP::THandle<CP::TFADCHit>(new CP::TFADCHit(hit)));

        double expectedCharge = 0.0;
        for (CP::THitSelection::iterator uHit = uHits.begin();
             uHit != uHits.end(); ++uHit) {
            expectedCharge += (*uHit)->GetCharge();
        }
        for (CP::THitSelection::iterator vHit = vHits.begin();
             vHit != vHits.end(); ++vHit) {
            expectedCharge += (*vHit)->GetCharge();
        }
        for (CP::THitSelection::iterator xHit = xHits.begin();
             xHit != xHits.end(); ++xHit) {
            expectedCharge += (*xHit)->GetCharge();
        }
        
        CP::THitSelection hits3D;
        for (CP::THitSelection::iterator uHit = uHits.begin();
             uHit != uHits.end(); ++uHit) {
            for (CP::THitSelection::iterator vHit = vHits.begin();
                 vHit != vHits.end(); ++vHit) {
                for (CP::THitSelection::iterator xHit = xHits.begin();
                     xHit != xHits.end(); ++xHit) {
                    CP::TWritableReconHit hit(*xHit,*vHit,*uHit);
                    hits3D.push_back(CP::THandle<CP::TReconHit>(
                                         new CP::TReconHit(hit)));
                }
            }
        }

        // Share the charge among the 3D hits so that the total charge in the
        // event is not overcounted.
        CP::TShareCharge share;

        // Fill the charge sharing object.
        for (CP::THitSelection::iterator h = hits3D.begin();
             h != hits3D.end(); ++h) {
            CP::ShareCharge::TMeasurementGroup& group = share.AddGroup(*h);
            CP::THandle<CP::TReconHit> groupHit = *h;
            for (int i=0; i<groupHit->GetConstituentCount(); ++i) {
                CP::THandle<CP::THit> hit = groupHit->GetConstituent(i); 
                group.AddMeasurement(hit, hit->GetCharge());
            }
        }

        share.Solve();

        double finalCharge = 0.0;
        for (CP::TShareCharge::Groups::const_iterator g
                 = share.GetGroups().begin();
             g != share.GetGroups().end(); ++g) {
            CP::THandle<CP::TWritableReconHit> groupHit = g->GetObject();
            int group = 0;
            finalCharge += g->GetGroupCharge();
            for(CP::ShareCharge::TLinks::const_iterator c
                    = g->GetLinks().begin();
                c != g->GetLinks().end(); ++c) {
                CP::THandle<CP::THit> hit = (*c)->GetMeasurement()->GetObject();
                CP::TGeometryId geomId = hit->GetGeomId();
                if (group == 0) {
                    group = CP::GeomId::Captain::GetWireNumber(geomId);
                }
                else if (group != CP::GeomId::Captain::GetWireNumber(geomId)) {
                    group = -1;
                }
            }
        }
        ensure_tolerance("Charge sharing perserves normalization.", 
                         finalCharge,
                         expectedCharge, 0.0001);
    }
};

// Local Variables:
// mode:c++
// c-basic-offset:4
// End:
