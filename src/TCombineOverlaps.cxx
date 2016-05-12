#include "TCombineOverlaps.hxx"
#include "CreateCluster.hxx"
#include "CreateShower.hxx"
#include "CreateTrack.hxx"
#include "TTrackFit.hxx"

#include <THandle.hxx>
#include <TReconTrack.hxx>
#include <TReconCluster.hxx>
#include <TCaptLog.hxx>
#include <TRuntimeParameters.hxx>

#include <memory>
#include <cmath>
#include <set>
#include <list>
#include <algorithm>

/// A set of 2D hits.  The type doesn't prevent 3D hits from being added,
/// but this set type should only hold 2D hits.
CP::TCombineOverlaps::TCombineOverlaps()
    : TAlgorithm("TCombineOverlaps", 
                 "Combine overlapping objects") {
    fOverlapCut = CP::TRuntimeParameters::Get().GetParameterD(
            "captRecon.combineOverlaps.overlap");
}

CP::TCombineOverlaps::~TCombineOverlaps() { }

double CP::TCombineOverlaps::CheckOverlap(
    CP::THandle<CP::TReconBase> object1,
    CP::THandle<CP::TReconBase> object2) {

    // Divide objects into 2D hits.
    std::set< CP::THandle<CP::THit> > set1u;
    std::set< CP::THandle<CP::THit> > set1v;
    std::set< CP::THandle<CP::THit> > set1x;
    for (CP::THitSelection::iterator h = object1->GetHits()->begin();
         h != object1->GetHits()->end(); ++h) {
        for (int i = 0; i < (*h)->GetConstituentCount(); ++i) {
            CP::THandle<CP::THit> w = (*h)->GetConstituent(i);
            CP::TGeometryId id = w->GetGeomId();
            if (CP::GeomId::Captain::IsUWire(id)) set1u.insert(w);
            if (CP::GeomId::Captain::IsVWire(id)) set1v.insert(w);
            if (CP::GeomId::Captain::IsXWire(id)) set1x.insert(w);
        }
    }

    std::set< CP::THandle<CP::THit> > set2u;
    std::set< CP::THandle<CP::THit> > set2v;
    std::set< CP::THandle<CP::THit> > set2x;
    for (CP::THitSelection::iterator h = object2->GetHits()->begin();
         h != object2->GetHits()->end(); ++h) {
        for (int i = 0; i < (*h)->GetConstituentCount(); ++i) {
            CP::THandle<CP::THit> w = (*h)->GetConstituent(i);
            CP::TGeometryId id = w->GetGeomId();
            if (CP::GeomId::Captain::IsUWire(id)) set2u.insert(w);
            if (CP::GeomId::Captain::IsVWire(id)) set2v.insert(w);
            if (CP::GeomId::Captain::IsXWire(id)) set2x.insert(w);
        }
    }


    // Subtract the larger from the smaller and find the overlap.
    double set1USize = set1u.size(); // YES, IT IS A DOUBLE.
    double set2USize = set2u.size(); // YES, IT IS A DOUBLE.
    std::set< CP::THandle<CP::THit> > uniqueU;
    double overlapU = 0.0;
    if (set1USize < set2USize) {
        std::set_difference(set1u.begin(), set1u.end(),
                            set2u.begin(), set2u.end(),
                            std::inserter(uniqueU, uniqueU.end()));
        overlapU = 1.0 - uniqueU.size()/set1USize;
    }
    else {
        std::set_difference(set2u.begin(), set2u.end(),
                            set1u.begin(), set1u.end(),
                            std::inserter(uniqueU, uniqueU.end()));
        overlapU = 1.0 - uniqueU.size()/set2USize;
    }

    CaptNamedInfo("Combine", "U Overlap: " << overlapU
                     << " out of " << set1USize 
                     << " (" << object1->GetUniqueID() << ")"
                     << " and " << set2USize
                     << " (" << object2->GetUniqueID() << ")");

    // Subtract the larger from the smaller and find the overlap.
    double set1VSize = set1v.size(); // YES, IT IS A DOUBLE.
    double set2VSize = set2v.size(); // YES, IT IS A DOUBLE.
    std::set< CP::THandle<CP::THit> > uniqueV;
    double overlapV = 0.0;
    if (set1VSize < set2VSize) {
        std::set_difference(set1v.begin(), set1v.end(),
                            set2v.begin(), set2v.end(),
                            std::inserter(uniqueV, uniqueV.end()));
        overlapV = 1.0 - uniqueV.size()/set1VSize;
    }
    else {
        std::set_difference(set2v.begin(), set2v.end(),
                            set1v.begin(), set1v.end(),
                            std::inserter(uniqueV, uniqueV.end()));
        overlapV = 1.0 - uniqueV.size()/set2VSize;
    }

    CaptNamedInfo("Combine", "V Overlap: " << overlapV
                  << " out of " << set1VSize 
                  << " (" << object1->GetUniqueID() << ")"
                  << " and " << set2VSize
                  << " (" << object2->GetUniqueID() << ")");
    
    // Subtract the larger from the smaller and find the overlap.
    double set1XSize = set1x.size(); // YES, IT IS A DOUBLE.
    double set2XSize = set2x.size(); // YES, IT IS A DOUBLE.
    std::set< CP::THandle<CP::THit> > uniqueX;
    double overlapX = 0.0;
    if (set1XSize < set2XSize) {
        std::set_difference(set1x.begin(), set1x.end(),
                            set2x.begin(), set2x.end(),
                            std::inserter(uniqueX, uniqueX.end()));
        overlapX = 1.0 - uniqueX.size()/set1XSize;
    }
    else {
        std::set_difference(set2x.begin(), set2x.end(),
                            set1x.begin(), set1x.end(),
                            std::inserter(uniqueX, uniqueX.end()));
        overlapX = 1.0 - uniqueX.size()/set2XSize;
    }

    CaptNamedInfo("Combine", "X Overlap: " << overlapX
                     << " out of " << set1XSize 
                     << " (" << object1->GetUniqueID() << ")"
                     << " and " << set2XSize
                     << " (" << object2->GetUniqueID() << ")");

    double overlap = overlapU;
    overlap = std::min(overlap,overlapV);
    overlap = std::min(overlap,overlapX);
    
    return overlap;
}

CP::THandle<CP::TReconBase> 
CP::TCombineOverlaps::MergeObjects(CP::THandle<CP::TReconBase> object1,
             CP::THandle<CP::TReconBase> object2) {
    std::set< CP::THandle<THit> > objectHits;

    // Insert the object hits into a set.  This will make sure they are unique
    // in the new object.
    CP::THitSelection::iterator begin = object1->GetHits()->begin();
    CP::THitSelection::iterator end = object1->GetHits()->end();
    while (begin != end) objectHits.insert(*(begin++));

    begin = object2->GetHits()->begin();
    end = object2->GetHits()->end();
    while (begin != end) objectHits.insert(*(begin++));

    CP::THandle<CP::TReconTrack> merged
        = CreateTrackFromHits("TCombineOverlaps",
                              objectHits.begin(), objectHits.end());
    
    if (!merged) {
        CaptError("No combined object constructed");
    }

    return merged;
}

CP::THandle<CP::TAlgorithmResult>
CP::TCombineOverlaps::Process(const CP::TAlgorithmResult& input,
                               const CP::TAlgorithmResult&,
                               const CP::TAlgorithmResult&) {
    
    CP::THandle<CP::TReconObjectContainer> inputObjects 
        = input.GetResultsContainer();

    CaptLog("TCombineOverlaps Process " << GetEvent().GetContext());

    if (!inputObjects) {
        CaptError("No input objects");
        return CP::THandle<CP::TAlgorithmResult>();
    }

    // Create the output containers.
    CP::THandle<CP::TAlgorithmResult> result = CreateResult();
    std::unique_ptr<CP::TReconObjectContainer> 
        final(new CP::TReconObjectContainer("final"));

    // Create a list to keep objects.  When the list is empty, all the objects
    // that need to be merge have been merge.  This is a pseudo stack where
    // objects are popped off of the top until it's empty, but they are also
    // removed from the center when a match is found.
    typedef std::list< CP::THandle<CP::TReconBase> > ObjectList;
    ObjectList objectList;

    // Make a copy of all the input objects to a list to be manipulated.
    std::copy(inputObjects->begin(), inputObjects->end(),
              std::back_inserter(objectList));

    // Pop an object off the stack and see if it should be merged.  If the track
    // is merged, the result is pushed back on the stack.  If it doesn't get
    // merge, the track get's pushed into the final object container.
    while (!objectList.empty()) {
        CP::THandle<CP::TReconBase> object1 = objectList.front();
        objectList.pop_front();
        CaptNamedInfo("Combine",
                      "Stack Size: " << objectList.size()
                      << "    Object size: " << object1->GetNodes().size()
                      << "    UID: " << object1->GetUniqueID());

        for (ObjectList::iterator t = objectList.begin();
             t!=objectList.end(); ++t) {
            CP::THandle<CP::TReconBase> object2 = *t;
            
            double match = CheckOverlap(object1,object2);

            CaptNamedVerbose("Combine", "Object overlap " << match 
                             << " w/ objects: " 
                             << object1->GetUniqueID()
                             << " ("  << object1->GetHits()->size() << " hits)"
                             << ", " << object2->GetUniqueID()
                             << "  (" << object2->GetHits()->size()
                             << " hits)");

            if (match < fOverlapCut) continue;
            
            ///////////////////////////////////////////////////////
            // If we get here then the tracks should be merged.
            ///////////////////////////////////////////////////////
            CP::THandle<CP::TReconBase> merged = MergeObjects(object1,object2);

            if (!merged) {
                CaptError("Objects not merged");
                continue;
            }

            // Remove the object iterator.
            objectList.erase(t);

            // Merge the tracks.
            CaptNamedInfo("Combine", "Matched: " << match << " -- "
                         << object1->GetUniqueID() 
                         << " w/ " << object1->GetHits()->size()
                         << ", " << object2->GetUniqueID()
                         << " w/ " << object2->GetHits()->size()
                         << " into " << merged->GetUniqueID()
                         << " w/ " << merged->GetHits()->size() << " hits");

            // Put the new track back into the list.
            objectList.push_front(merged);

            // Clear out the track variable.
            object1 = CP::THandle<CP::TReconBase>();

            // Don't continue the loop looking for a pair of objects to
            // merge.
            break;
        }

        // If we get to the bottom of the loop looking for a pair of objects
        // to merge then there is nothing to merge with this object.  Push it
        // on to the final objects.
        if (object1) {
            CaptNamedInfo("Combine", "Save Object UID: "
                          << object1->GetUniqueID());

            // Check if this is a track and it should be fit.
            CP::THandle<CP::TReconTrack> track = object1;
            if (track) {
                if (!track->CheckStatus(CP::TReconBase::kSuccess)) {
                    CaptNamedInfo("Combine", "Fit track UID "
                                  << track->GetUniqueID());
                    TTrackFit fitter;
                    object1 = fitter(track);
                }
            }

            final->push_back(object1);
        }

    }

    result->AddResultsContainer(final.release());

    return result;
}
