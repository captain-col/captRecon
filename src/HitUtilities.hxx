#ifndef HitUtilities_hxx_seen

#include <THit.hxx>
#include <THandle.hxx>
#include <THitSelection.hxx>
#include <TReconBase.hxx>

namespace CP {
namespace hits {
    /// Return true if two hit handles are equal.  Hits are equal if they
    /// point to the same address, and is equivalent to the
    /// CP::THandle<CP::THit>::operator==().  They are also equal if they have
    /// constituents, and all constituents are shared.
    bool Equal(const CP::THandle<CP::THit>& a, const CP::THandle<CP::THit>& b);

    /// Remove hits from the first hit selection that have an equal entry
    /// in the second hit selection.
    void Subtract(CP::THitSelection& a, const CP::THitSelection& b);
    
    /// Remove duplicate hits from a hit selection.
    void Unique(CP::THitSelection& a);

    /// Collect all of the hits used by TReconBase objects in a
    /// reconstruction object container into a single hit selection.
    template<typename T>
    CP::THandle<CP::THitSelection> 
    ReconHits(T begin, T end) {
        CP::THandle<CP::THitSelection> hits(new CP::THitSelection);
        
        while (begin != end) {
            // Add the hits for the object.
            CP::THandle<CP::THitSelection> objHits = (*begin)->GetHits();
            if (objHits) {
                for (CP::THitSelection::iterator hit = objHits->begin();
                     hit != objHits->end();
                     ++hit) {
                    hits->AddHit(*hit);
                }
            }
            // Add hits for the constituents.  
            CP::THandle<CP::TReconObjectContainer> parts 
                = (*begin)->GetConstituents();
            if (parts) {
                objHits = ReconHits(parts->begin(), parts->end());
                if (objHits) {
                    for (CP::THitSelection::iterator hit = objHits->begin();
                         hit != objHits->end();
                         ++hit) {
                        hits->AddHit(*hit);
                    }
                }
            }
            ++begin;
        }
        
        return hits;
    }

};
};
#endif
