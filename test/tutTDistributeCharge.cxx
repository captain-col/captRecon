#include <TDistributeCharge.hxx>

#include <HEPUnits.hxx>
#include <TCaptLog.hxx>
#include <TFADCHit.hxx>
#include <TReconHit.hxx>
#include <CaptGeomId.hxx>

#include <tut.h>

namespace tut {
    struct baseDistributeCharge {
        baseDistributeCharge() {
            // Run before each test.
        }
        ~baseDistributeCharge() {
            // Run after each test.
        }
    };

    // Declare the test
    typedef test_group<baseDistributeCharge>::object testDistributeCharge;
    test_group<baseDistributeCharge> groupDistributeCharge("TDistributeCharge");

    // Test the declaration.
    template<> template<> void testDistributeCharge::test<1> () {
        CP::TDistributeCharge distribute;
    }

    // Test the declaration.
    template<> template<> void testDistributeCharge::test<2> () {
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
        for (CP::THitSelection::iterator h = xHits.begin();
             h != xHits.end(); ++h) {
            expectedCharge += (*h)->GetCharge();
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

        // Distribute the charge among the 3D hits so that the total charge in
        // the event is not overcounted.
        CP::TDistributeCharge distribute;

        // Fill the charge sharing object.
        for (CP::THitSelection::iterator h = hits3D.begin();
             h != hits3D.end(); ++h) {
            CP::DistributeCharge::TMeasurementGroup& group
                = distribute.AddGroup(*h);
            CP::THandle<CP::TReconHit> groupHit = *h;
            CP::THandle<CP::THit> xhit = groupHit->GetConstituent(0);
            CP::THandle<CP::THit> vhit = groupHit->GetConstituent(1);
            CP::THandle<CP::THit> uhit = groupHit->GetConstituent(2);
            double physicsWeight = 1.0;
            for (int i=0; i<groupHit->GetConstituentCount(); ++i) {
                CP::THandle<CP::THit> hit = groupHit->GetConstituent(i); 
                group.AddMeasurement(hit, hit->GetCharge(), physicsWeight);
            }
        }

        distribute.Solve();

        double finalCharge = 0.0;
        for (CP::TDistributeCharge::Groups::const_iterator g
                 = distribute.GetGroups().begin();
             g != distribute.GetGroups().end(); ++g) {
            CP::THandle<CP::TWritableReconHit> groupHit = g->GetObject();
            double gc =  g->GetGroupCharge();
            finalCharge += gc;
#ifdef DEBUG_OUTPUT_ENABLED
#undef DEBUG_OUTPUT_ENABLED
            std::cout << "gc " << gc;
            for(CP::DistributeCharge::TLinks::const_iterator c
                    = g->GetLinks().begin();
                c != g->GetLinks().end(); ++c) {
                CP::THandle<CP::THit> hit = (*c)->GetMeasurement()->GetObject();
                CP::TGeometryId geomId = hit->GetGeomId();
                std::cout << " (" << CP::GeomId::Captain::GetWireNumber(geomId)
                          << "," << hit->GetCharge()
                          << "," << (*c)->GetPhysicsWeight()
                          << "," << (*c)->GetCharge() << ")";
            }
            std::cout << std::endl;
#endif
        }

        ensure_tolerance("Charge distribution perserves normalization.", 
                         finalCharge,
                         expectedCharge, 0.0001);
    }

    // Test the declaration.
    template<> template<> void testDistributeCharge::test<3> () {
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
        hit.SetCharge(1.0);
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
        hit.SetCharge(1.0);
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
        hit.SetCharge(5.0);
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
        for (CP::THitSelection::iterator h = xHits.begin();
             h != xHits.end(); ++h) {
            expectedCharge += (*h)->GetCharge();
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

        // Distribute the charge among the 3D hits so that the total charge in
        // the event is not overcounted.
        CP::TDistributeCharge distribute;

        // Fill the charge sharing object.
        for (CP::THitSelection::iterator h = hits3D.begin();
             h != hits3D.end(); ++h) {
            CP::DistributeCharge::TMeasurementGroup& group
                = distribute.AddGroup(*h);
            CP::THandle<CP::TReconHit> groupHit = *h;
            double physicsWeight = 1.0;
            for (int i=0; i<groupHit->GetConstituentCount(); ++i) {
                CP::THandle<CP::THit> hit = groupHit->GetConstituent(i); 
                group.AddMeasurement(hit, hit->GetCharge(), physicsWeight);
            }
        }

        distribute.Solve();

        double finalCharge = 0.0;
        for (CP::TDistributeCharge::Groups::const_iterator g
                 = distribute.GetGroups().begin();
             g != distribute.GetGroups().end(); ++g) {
            CP::THandle<CP::TWritableReconHit> groupHit = g->GetObject();
            double gc =  g->GetGroupCharge();
            finalCharge += gc;
#ifdef DEBUG_OUTPUT_ENABLED
#undef DEBUG_OUTPUT_ENABLED
            std::cout << "gc " << gc;
            for(CP::DistributeCharge::TLinks::const_iterator c
                    = g->GetLinks().begin();
                c != g->GetLinks().end(); ++c) {
                CP::THandle<CP::THit> hit = (*c)->GetMeasurement()->GetObject();
                CP::TGeometryId geomId = hit->GetGeomId();
                std::cout << " (" << CP::GeomId::Captain::GetWireNumber(geomId)
                          << "," << hit->GetCharge()
                          << "," << (*c)->GetPhysicsWeight()
                          << "," << (*c)->GetCharge() << ")";
            }
            std::cout << std::endl;
#endif
        }

        ensure_tolerance("Charge distribution perserves normalization.", 
                         finalCharge,
                         expectedCharge, 0.0001);
    }

    // Test the declaration.
    template<> template<> void testDistributeCharge::test<4> () {
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
        hit.SetCharge(1.0);
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
        hit.SetCharge(1.0);
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
        hit.SetCharge(5.0);
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
        for (CP::THitSelection::iterator h = xHits.begin();
             h != xHits.end(); ++h) {
            expectedCharge += (*h)->GetCharge();
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

        // Distribute the charge among the 3D hits so that the total charge in
        // the event is not overcounted.
        CP::TDistributeCharge distribute;

        // Fill the charge sharing object.
        for (CP::THitSelection::iterator h = hits3D.begin();
             h != hits3D.end(); ++h) {
            CP::DistributeCharge::TMeasurementGroup& group
                = distribute.AddGroup(*h);
            CP::THandle<CP::TReconHit> groupHit = *h;
            CP::THandle<CP::THit> xhit = groupHit->GetConstituent(0);
            CP::THandle<CP::THit> vhit = groupHit->GetConstituent(1);
            CP::THandle<CP::THit> uhit = groupHit->GetConstituent(2);
            double physicsWeight = 1.0;
            if (CP::GeomId::Captain::GetWireNumber(xhit->GetGeomId())
                != CP::GeomId::Captain::GetWireNumber(vhit->GetGeomId())) {
                physicsWeight *= 0.1;
            }
            if (CP::GeomId::Captain::GetWireNumber(xhit->GetGeomId())
                != CP::GeomId::Captain::GetWireNumber(uhit->GetGeomId())) {
                physicsWeight *= 0.1;
            }
            if (CP::GeomId::Captain::GetWireNumber(uhit->GetGeomId())
                != CP::GeomId::Captain::GetWireNumber(vhit->GetGeomId())) {
                physicsWeight *= 0.1;
            }
            for (int i=0; i<groupHit->GetConstituentCount(); ++i) {
                CP::THandle<CP::THit> hit = groupHit->GetConstituent(i); 
                group.AddMeasurement(hit, hit->GetCharge(), physicsWeight);
            }
        }

        distribute.Solve();

        double finalCharge = 0.0;
        for (CP::TDistributeCharge::Groups::const_iterator g
                 = distribute.GetGroups().begin();
             g != distribute.GetGroups().end(); ++g) {
            CP::THandle<CP::TWritableReconHit> groupHit = g->GetObject();
            double gc =  g->GetGroupCharge();
            finalCharge += gc;
#ifdef DEBUG_OUTPUT_ENABLED
#undef DEBUG_OUTPUT_ENABLED
            std::cout << "gc " << gc;
            for(CP::DistributeCharge::TLinks::const_iterator c
                    = g->GetLinks().begin();
                c != g->GetLinks().end(); ++c) {
                CP::THandle<CP::THit> hit = (*c)->GetMeasurement()->GetObject();
                CP::TGeometryId geomId = hit->GetGeomId();
                std::cout << " (" << CP::GeomId::Captain::GetWireNumber(geomId)
                          << "," << hit->GetCharge()
                          << "," << (*c)->GetPhysicsWeight()
                          << "," << (*c)->GetCharge() << ")";
            }
            std::cout << std::endl;
#endif
        }

        ensure_tolerance("Charge distribution perserves normalization.", 
                         finalCharge,
                         expectedCharge, 0.0001);
    }
};

// Local Variables:
// mode:c++
// c-basic-offset:4
// End:
