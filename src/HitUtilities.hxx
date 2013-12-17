#ifndef HitUtilities_hxx_seen
#define HitUtilities_hxx_seen

#include <THit.hxx>
#include <THandle.hxx>
#include <THitSelection.hxx>
#include <TReconBase.hxx>
#include <TReconNode.hxx>

#include <set>

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

    /// Collect all of the hits used in a TReconBase object into a single set.
    void ReconHits(CP::THandle< CP::TReconBase > object, 
                    std::set< CP::THandle<CP::THit> >& output);

    /// Collect all of the hits used by TReconBase objects in a
    /// reconstruction object container into a single set of hits.
    template<typename InputIterator>
    void ReconHits(InputIterator begin, InputIterator end, 
                   std::set< CP::THandle<CP::THit> >& output) {
        while (begin != end) {
            ReconHits(*begin, output);
            ++begin;
        }
    }

    /// Collect all of the hits used by TReconBase objects in a reconstruction
    /// object container into a THitSelection.  A hit will only appear in the
    /// selection once.
    template<typename InputIterator>
    CP::THandle< CP::THitSelection >
    ReconHits(InputIterator begin, InputIterator end) {
        std::set< CP::THandle<CP::THit> > theSet;
        ReconHits(begin,end,theSet);

        CP::THandle<CP::THitSelection> hits(new CP::THitSelection);
        hits->reserve(theSet.size());
        std::copy(theSet.begin(), theSet.end(), std::back_inserter(*hits));

        return hits;
    }

};
};
#endif
