#include "TRemoveOutliers.hxx"
#include "TCaptLog.hxx"

#include <CaptGeomId.hxx>

#include <set>
#include <cmath>

CP::TRemoveOutliers::TRemoveOutliers(int maxSize, double small, double big) 
    : fMaxOutlierSize(maxSize), 
      fSmallGroupFraction(small), 
      fBigGroupFraction(big) {} 

bool CP::TRemoveOutliers::IsOutlier(CP::THandle<CP::THit> hit, int index) {
    CP::THandle<CP::TReconHit> reconHit = hit;
    HitInfo& info = fHitMap[reconHit->GetConstituent(index)];

    int plane = CP::GeomId::Captain::GetWirePlane(
        reconHit->GetConstituent(index)->GetGeomId());
    // The number of the first crossing wire in the hit.
    int number1 = CP::GeomId::Captain::GetWireNumber(
        reconHit->GetConstituent((index+1)%3)->GetGeomId());
    // The number of the second crossing wire in the hit.
    int number2 = CP::GeomId::Captain::GetWireNumber(
        reconHit->GetConstituent((index+2)%3)->GetGeomId());

    // Sets of the crossing wires.  These correspond to the planes for number1
    // and number2.
    std::set<int> set1;
    std::set<int> set2;
    for (CP::THitSelection::iterator h = info.fContainedBy.begin();
         h != info.fContainedBy.end(); ++h) {
        CP::THandle<CP::TReconHit> neighbor = *h;
        for (int i=0; i<neighbor->GetConstituentCount(); ++i) {
            CP::TGeometryId geomId 
                = neighbor->GetConstituent(i)->GetGeomId();
            int np = CP::GeomId::Captain::GetWirePlane(geomId);
            if (plane == np) continue;
            // Wire set will be 0 or 1 and correspondes to the other two
            // planes that are NOT the plane for the current input hit.
            int wireSet = (2+np-plane) % 3;
            int wire = CP::GeomId::Captain::GetWireNumber(geomId);
            if (!wireSet) set1.insert(wire);
            else set2.insert(wire);
        }
    }

    double delta1 = LocalGroup(number1, set1);
    double delta2 = LocalGroup(number2, set2);

    if (delta1 > fMaxOutlierSize && delta2 > fMaxOutlierSize) return false;

    delta1 /= set1.size();
    delta2 /= set2.size();

    if (delta1 < fSmallGroupFraction) return true;
    if (delta2 < fSmallGroupFraction) return true;

    if (delta1 < fBigGroupFraction && delta2 < fBigGroupFraction) return true;

    return false;
}

int CP::TRemoveOutliers::LocalGroup(int wire, std::set<int>& input) {
    // Find the current hit wire in the set.
    std::set<int>::iterator current = input.find(wire);
    std::set<int>::iterator first = current;
    std::set<int>::iterator next = first;
    int count = 1;
    while (next != input.begin()) {
        --next;
        if (std::abs(*first-*next) > 3) break;
        first = next;
        ++count;
    }
    std::set<int>::iterator last = current;
    next = last;
    while (true) {
        ++next;
        if (next == input.end()) break;
        if (std::abs(*last-*next) > 2) break;
        last = next;
        ++count;
    }

    return count;
}

void CP::TRemoveOutliers::Apply(CP::THitSelection& hits) {
    fHitMap.clear();

    for (CP::THitSelection::iterator h = hits.begin();
         h != hits.end(); ++h) {
        CP::THandle<CP::TReconHit> reconHit = *h;
        if (!reconHit) {
            CaptError("Input hit is not a TReconHit");
            throw;
        }
        // Fill the one to many map of 2D Hits to 3D hits. 
        for (int i=0; i<reconHit->GetConstituentCount(); ++i) {
            CP::THandle<CP::THit> hit = reconHit->GetConstituent(i);
            fHitMap[hit].fContainedBy.push_back(reconHit);
        }
    }
 
    CP::THitSelection::iterator last = hits.begin();
    for (CP::THitSelection::iterator h = hits.begin();
         h != hits.end(); ++h) {
        CP::THandle<CP::TReconHit> reconHit = *h;
        int outlierCount = 0;
        for (int i=0; i<reconHit->GetConstituentCount(); ++i) {
            if (IsOutlier(reconHit,i)) ++outlierCount;
        }
        if (outlierCount>1) continue;
        *last = *h;
        ++last;
    }

    hits.erase(last,hits.end());
}
