#include <algorithm>
#include <iterator>

#include "HitUtilities.hxx"

#include <THitSelection.hxx>
#include <THandle.hxx>
#include <HEPUnits.hxx>

namespace {
    // Take a hit and unpack it into it's most basic constituents.  For a
    // "normal" CP::THit, this just puts the hit into the result.  For a hit
    // with constituents, this unpacks each individual constiuent and puts the
    // the low level hits into the result.
    template <class T>
    void UnpackHits(const CP::THandle<CP::THit>& h, T result) {
        int cCount = h->GetConstituentCount();
        if (cCount < 1) {
            *result++ = h;
            return;
        }
        for (int i = 0; i<cCount; ++i) {
            UnpackHits(h->GetConstituent(i),result);
        }
    }
};

bool CP::hits::Equal(const CP::THandle<CP::THit>& a, 
                     const CP::THandle<CP::THit>& b) {
    // Check if the two handles are pointing to the same object.
    if (CP::GetPointer(a) == CP::GetPointer(b)) return true;
#ifndef DEEP_EQUALITY_CHECK
    return false;
#else
    // The handles point to different objects, but they might still be the
    // same "hit".  The hits will be equal if they have the same constituents. 
    if (a->GetConstituentCount() != b->GetConstituentCount()) return false;

    // Neither hit has constituents.  The two hits will be equal if they are
    // for the same object, and have "identical" times.
    if (a->GetConstituentCount() == 0 && b->GetConstituentCount() == 0) {
        if (a->GetGeomId() != b->GetGeomId()) return false;
        double dt = std::abs(a->GetTime() - b->GetTime());
        if (dt > 1*unit::nanosecond) return false;
        return true;
    }

    // Unpack the hits into the most basic hits.
    CP::THitSelection aHits;
    UnpackHits(a,std::back_inserter(aHits));
    
    CP::THitSelection bHits;
    UnpackHits(b,std::back_inserter(bHits));

    // Check that none of the constituent hits are not equal.
    for (CP::THitSelection::iterator m=aHits.begin(); m!=aHits.end(); ++m) {
        for (CP::THitSelection::iterator n=bHits.begin(); n!=bHits.end(); ++n) {
            if (!Equal(*m,*n)) return false;
        }
    }

    return true;
#endif
}
    
void CP::hits::Subtract(CP::THitSelection& a, const CP::THitSelection& b) {
    CP::THitSelection::const_iterator begin = a.begin();
    CP::THitSelection::const_iterator end = a.end();
    CP::THitSelection::iterator result = a.begin();
    while (begin != end) {
        bool found = false;
        for (CP::THitSelection::const_iterator h = b.begin();
             h != b.end(); ++h) {
            if (Equal(*begin,*h)) {
                found = true;
                break;
            }
        }
        if (!found) {
            CP::THandle<CP::THit> copy(*begin);
            *result = copy;
            ++result;
        }
        ++begin;
    }
    if (result!=a.end()) a.erase(result,a.end());
}

void CP::hits::Unique(CP::THitSelection& a) {
    CP::THitSelection::const_iterator begin = a.begin();
    CP::THitSelection::const_iterator end = a.end();
    CP::THitSelection::iterator result = a.begin();
    while (begin != end) {
        bool found = false;
        for (CP::THitSelection::const_iterator h = a.begin();
             h != begin; ++h) {
            if (Equal(*begin,*h)) {
                found = true;
                break;
            }
        }
        if (!found) {
            CP::THandle<CP::THit> copy(*begin);
            *result = copy;
            ++result;
        }
        ++begin;
    }
    if (result!=a.end()) a.erase(result,a.end());
}

CP::THandle<CP::THitSelection> 
CP::hits::ReconHits(const CP::TReconObjectContainer& input) {
    CP::THandle<CP::THitSelection> hits(new CP::THitSelection);

    for (CP::TReconObjectContainer::const_iterator o = input.begin();
         o != input.end();
         ++o) {
        // Add the hits for the object.
        CP::THandle<CP::THitSelection> objHits = (*o)->GetHits();
        if (objHits) {
            for (CP::THitSelection::iterator hit = objHits->begin();
                 hit != objHits->end();
                 ++hit) {
                hits->AddHit(*hit);
            }
        }
        // Add hits for the constituents.  
        CP::THandle<CP::TReconObjectContainer> parts = (*o)->GetConstituents();
        if (parts) {
            objHits = ReconHits(*parts);
            if (objHits) {
                for (CP::THitSelection::iterator hit = objHits->begin();
                     hit != objHits->end();
                     ++hit) {
                    hits->AddHit(*hit);
                }
            }
        }
    }

    return hits;
}
