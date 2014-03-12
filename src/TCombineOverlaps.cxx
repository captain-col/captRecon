#include "TCombineOverlaps.hxx"
#include "CreateCluster.hxx"
#include "CreateShower.hxx"

#include <THandle.hxx>
#include <TReconTrack.hxx>
#include <TReconCluster.hxx>
#include <TCaptLog.hxx>
 
#include <memory>
#include <cmath>
#include <set>
#include <list>
#include <algorithm>

/// A set of 2D hits.  The type doesn't prevent 3D hits from being added,
/// but this set type should only hold 2D hits.
typedef std::set< CP::THandle< CP::THit > > HitSet;

/// A structrue to hold the 2D hits that contribute to a reconstruction
/// object.
typedef struct {
    HitSet fXHits;
    HitSet fVHits;
    HitSet fUHits;
} XVUHits;


CP::TCombineOverlaps::TCombineOverlaps()
    : TAlgorithm("TCombineOverlaps", 
                 "Combine overlapping objects") {
    fOverlapCut = 0.9;
}

CP::TCombineOverlaps::~TCombineOverlaps() { }

double CP::TCombineOverlaps::CheckOverlap(
    CP::THandle<CP::TReconBase> object1,
    CP::THandle<CP::TReconBase> object2) {

    // Divide into objects 2D hits.
    std::set< CP::THandle<CP::THit> > set1;
    for (CP::THitSelection::iterator h = object1->GetHits()->begin();
         h != object1->GetHits()->end(); ++h) {
        for (int i = 0; i < (*h)->GetConstituentCount(); ++i) {
            CP::THandle<CP::THit> w = (*h)->GetConstituent(i);
            set1.insert(w);
        }
    }
    double set1Size = set1.size(); // YES, IT IS A DOUBLE

    std::set< CP::THandle<CP::THit> > set2;
    for (CP::THitSelection::iterator h = object2->GetHits()->begin();
         h != object2->GetHits()->end(); ++h) {
        for (int i = 0; i < (*h)->GetConstituentCount(); ++i) {
            CP::THandle<CP::THit> w = (*h)->GetConstituent(i);
            set2.insert(w);
        }
    }
    double set2Size = set2.size(); // YES, IT IS A DOUBLE

    // Subtract the larger from the smaller.
    std::set< CP::THandle<CP::THit> > unique;
    double overlap;
    if (set1Size < set2Size) {
        std::set_difference(set1.begin(), set1.end(),
                            set2.begin(), set2.end(),
                            std::inserter(unique, unique.end()));
        overlap = 1.0 - unique.size()/set1Size;
    }
    else {
        std::set_difference(set2.begin(), set2.end(),
                            set1.begin(), set1.end(),
                            std::inserter(unique, unique.end()));
        overlap = 1.0 - unique.size()/set2Size;
    }

    CaptNamedVerbose("Combine", "Overlap: " << overlap 
                     << " out of " << set1Size 
                     << " (" << object1->GetUniqueID() << ")"
                     << " and " << set2Size
                     << " (" << object2->GetUniqueID() << ")");

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

    // Create a shower from the combined hits.
    CP::THandle<CP::TReconCluster> cluster 
        = CreateCluster("cluster", objectHits.begin(), objectHits.end());

    if (!cluster) {
        CaptError("No cluster constructed: " << objectHits.size());
        return cluster;
    }

    CP::THandle<CP::TReconShower> merged
        = CreateShower("shower", objectHits.begin(), objectHits.end(),
                       cluster->GetLongAxis());

    if (!merged) {
        CaptError("No shower constructed");
        return cluster;
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
    std::auto_ptr<CP::TReconObjectContainer> 
        final(new CP::TReconObjectContainer("final"));

    // Create a list to keep objects.  When the list is empty, all the objects
    // that need to be merge have been merge.  This is a pseudo stack where
    // objects are popped off of the top until it's empty, but they are also
    // removed from the center when a match is found.
    typedef std::list< CP::THandle<CP::TReconBase> > ObjectList;
    ObjectList objectList;

    // Make a copy of all of the tracks in the input objects, and save the
    // non-tracks to the output.
    std::copy(inputObjects->begin(), inputObjects->end(),
              std::back_inserter(objectList));


    // Pop a track off the stack and see if it should be merged.  If the track
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
            
            CaptNamedVerbose("Combine", "Check Objects"
                          << " w/ stack: " << objectList.size() 
                          << "  objects: " << object1->GetUniqueID()
                          << ", " << object2->GetUniqueID()
                          << "  sizes: " 
                          << object1->GetHits()->size() 
                          << ", " << object2->GetHits()->size());

            double match = CheckOverlap(object1,object2);
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
        // to merge with an unmerged object, then push it on to the final
        // objects.
        if (object1) {
            CaptNamedInfo("Combine", "Save Object");
            final->push_back(object1);
        }

    }

    result->AddResultsContainer(final.release());

    return result;
}
