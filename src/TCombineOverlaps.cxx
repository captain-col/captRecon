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

double CP::TCombineOverlaps::CountSetOverlaps(
    const CP::TCombineOverlaps::HitSet& set1,
    const CP::TCombineOverlaps::HitSet& set2) const {
    double set1Size = set1.size();
    double set2Size = set2.size();
    if (set1Size<1) return 1;
    if (set2Size<1) return 1;
    double overlap = 0;
    if (set1Size < set2Size) {
        for (HitSet::const_iterator h = set1.begin(); h!=set1.end(); ++h) {
            if (set2.find(*h) != set2.end()) overlap += 1.0;
        }
        return overlap/set1Size;
    }
    else {
        for (HitSet::const_iterator h = set2.begin(); h!=set2.end(); ++h) {
            if (set1.find(*h) != set1.end()) overlap += 1.0;
        }
        return overlap/set2Size;
    }
    return 0;
}
                                          
CP::THandle<CP::TReconBase> 
CP::TCombineOverlaps::MergeObjects(CP::THandle<CP::TReconBase> object1,
             CP::THandle<CP::TReconBase> object2) const { 
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

        // Divide objects into 2D hits.
        CP::TCombineOverlaps::HitSet set1u;
        CP::TCombineOverlaps::HitSet set1v;
        CP::TCombineOverlaps::HitSet set1x;
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

        for (ObjectList::iterator t = objectList.begin();
             t!=objectList.end(); ++t) {
            CP::THandle<CP::TReconBase> object2 = *t;

            CP::TCombineOverlaps::HitSet set2u;
            CP::TCombineOverlaps::HitSet set2v;
            CP::TCombineOverlaps::HitSet set2x;
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
            double overlapU = CountSetOverlaps(set1u,set2u);
            CaptNamedVerbose("Combine", "U Overlap: " << overlapU
                          << " (" << object1->GetUniqueID() << ")"
                          << " (" << object2->GetUniqueID() << ")");
            
            // Subtract the larger from the smaller and find the overlap.
            double overlapV = CountSetOverlaps(set1v,set2v);
            CaptNamedVerbose("Combine", "V Overlap: " << overlapV
                          << " (" << object1->GetUniqueID() << ")"
                          << " (" << object2->GetUniqueID() << ")");
            
            // Subtract the larger from the smaller and find the overlap.
            double overlapX = CountSetOverlaps(set1x,set2x);
            CaptNamedVerbose("Combine", "X Overlap: " << overlapX
                          << " (" << object1->GetUniqueID() << ")"
                          << " (" << object2->GetUniqueID() << ")");
            
            double overlap = overlapU;
            overlap = std::min(overlap,overlapV);
            overlap = std::min(overlap,overlapX);
            
            CaptNamedInfo("Combine", "Object overlap " << overlap
                             << " w/ objects: " 
                             << object1->GetUniqueID()
                             << " ("  << object1->GetHits()->size() << " hits)"
                             << ", " << object2->GetUniqueID()
                             << "  (" << object2->GetHits()->size()
                             << " hits)");

            if (overlap < fOverlapCut) continue;
            
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
            CaptNamedInfo("Combine", "Matched: " << overlap << " -- "
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

            // Don't continue the loop looking for a pair of objects to merge.
            // If the objects aren't to be merged, then the "continue" above
            // keeps the loop going.  If the objects were merged then the
            // result is back on the stack and we need to start the loop
            // again.
            break;
        }

        // If we get to the bottom of the loop looking for a pair of objects
        // to merge with a valid "object1, then there is nothing to be merged
        // with this object.  Push it on to the final objects.
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
