#ifndef CreateClusters_hxx_seen
#define CreateClusters_hxx_seen

#include "CreateCluster.hxx"

#include <ECore.hxx>
#include <TReconCluster.hxx>
#include <TReconShower.hxx>
#include <THandle.hxx>
#include <THit.hxx>
#include <HEPUnits.hxx>
#include <TCaptLog.hxx>

#include <algorithm>
#include <vector>
#include <list>

namespace CP {

    /// A base exception for the create track template.
    EXCEPTION(ECreateClusters,ECaptRecon);

    /// Create a CP::THandle<TReconObjectContainer> that contains clusters
    /// created from the input hits.  This needs the iterators to the range of
    /// hits to cluster, as well as an appropriate direction along which to
    /// cluster the hits.
    template<typename hitIterator>
    CP::THandle<CP::TReconObjectContainer> 
    CreateClusters(const char* name, hitIterator begin, hitIterator end,
                   const TVector3& approxDir,
                   double minLength, double maxLength,
                   double minSize, double maxSize);

    /// Create a CP::THandle<TReconObjectContainer> that contains clusters
    /// created from the input hits optimized for a shower.
    template<typename hitIterator>
    CP::THandle<CP::TReconObjectContainer> 
    CreateShowerClusters(const char* name, hitIterator begin, hitIterator end,
                         const TVector3& approxDir);

    /// Create a CP::THandle<TReconObjectContainer> that contains clusters
    /// created from the input hits optimized for a track.
    template<typename hitIterator>
    CP::THandle<CP::TReconObjectContainer> 
    CreateTrackClusters(const char* name, hitIterator begin, hitIterator end,
                        const TVector3& approxDir);

};

//////////////////////////////////////////////////////////////////
// IMPLEMENTATION
//////////////////////////////////////////////////////////////////

namespace {
    struct CreateClustersDirCompare {
        explicit CreateClustersDirCompare(const TVector3& dir) 
            : fDir(dir) {}
        bool operator() (const CP::THandle<CP::THit>& lhs,
                         const CP::THandle<CP::THit>& rhs) {
            return fDir*lhs->GetPosition() < fDir*rhs->GetPosition();
            
        }
    private:
        TVector3 fDir;
    };

    struct CreateClustersDirDiff {
        explicit CreateClustersDirDiff(const TVector3& dir): fDir(dir){}
        double operator() (const CP::THandle<CP::THit>& lhs,
                           const CP::THandle<CP::THit>& rhs) { 
#define DEBUG_CreateClustersDirDiff
#ifdef DEBUG_CreateClustersDirDiff
            if (!lhs) {
                CaptError("Invalid lhs");
                throw CP::ECreateClusters();
            }
            if (!rhs) {
                CaptError("Invalid rhs");
                throw CP::ECreateClusters();
            }
#endif
            return fDir*(lhs->GetPosition() - rhs->GetPosition());
        }
    private:
        TVector3 fDir;
    };
};
    
template<typename hitIterator>
CP::THandle<CP::TReconObjectContainer> 
CP::CreateClusters(const char* name, hitIterator begin, hitIterator end,
                   const TVector3& approxDir,
                   const double minLength, const double maxLength, 
                   const int minSize, const int maxSize) {
    
    TVector3 dir = approxDir;
    dir = dir.Unit();

    // Make the output object.
    CP::THandle<CP::TReconObjectContainer> 
        output(new CP::TReconObjectContainer(name));

    // copy the hits into local storage so they can be sorted.
    CP::THitSelection hits;
    hits.reserve(end-begin);
    std::copy(begin, end, std::back_inserter(hits));
    std::sort(hits.begin(), hits.end());
    CP::THitSelection::iterator kill = std::unique(hits.begin(), hits.end());
    hits.erase(kill,hits.end());

    // Sort the hits along the direction of the shower.
    std::sort(hits.begin(), hits.end(), CreateClustersDirCompare(dir));

    // Set the maximum length of a gap along the approxDir in a cluster before
    // the hits are enforced to be in separate clusters.
    const double maxGap = 10*unit::mm;

    CreateClustersDirDiff dirDiff(dir);
    
    // Record all of the places that the input selection should be split.  The
    // iterators should point to the first hit of each new cluster.  This
    // automatically sets the begin and end.
    typedef std::list<CP::THitSelection::iterator> Splits;
    Splits splits;
    splits.push_back(hits.begin());
    splits.push_back(hits.end());

    // Iterators to work with.
    CP::THitSelection::iterator curr = hits.begin();
    CP::THitSelection::iterator last = curr;
    CP::THitSelection::iterator first = curr;

    // Check to see if there are any gaps where the hits should be split.
    // This makes sure that if there is a gap in the THitSelection (along the
    // main direction), the gap isn't contained in a cluster.
    Splits::iterator split = splits.begin();
    while (curr != hits.end()) {
        first = curr;
        ++curr;
        // Make sure there are enough hits at the end to allow a split.
        if (hits.end()-curr < minSize) break;
        // Check there are enough hits since the last split to allow a new one.
        if (curr-last < minSize) continue;
        // Check if there needs to be a split.
        if (dirDiff(*curr, *first) < maxGap) continue;
        Splits::iterator next = split;
        split = splits.insert(++next, curr);
        last = curr;
    }

    Splits::iterator sBegin = splits.begin();
    Splits::iterator sEnd=sBegin; ++sEnd;
    while (sEnd != splits.end()) {
        // Not enough hits for a split.
        if ((*sEnd)-(*sBegin) < minSize) {
            sEnd = ++sBegin; ++sEnd;
            continue;
        }

        // Find the length of the current segment.
        double deltaS = dirDiff(*((*sEnd)-1), **sBegin);
        if (deltaS <= minLength) {
            sEnd = ++sBegin; ++sEnd;
            continue;
        }
        
        // Find the target length to evenly divide the region, while still
        // being less than maxLength.
        int nSegments = deltaS/maxLength + 1;
        double targetLength = deltaS/nSegments;
        int targetSize = (*sEnd - *sBegin)/maxSize + 1;
        
        // Check to see if there needs to be a split for the current segment.
        curr = *sBegin;
        first = curr;
        while (curr != *sEnd) {
            
            // Make sure each cluster has at least clustSize hits.
            if (curr-first < minSize) {
                ++curr;
                continue;
            }

            // Make sure the cluster isn't too long.
            double deltaS = dirDiff(*curr, *first);
            
            // Make sure each cluster is alt least minLength long.
            if (deltaS < minLength) {
                ++curr;
                continue;
            }

            // Make sure that the cluster at least the expected cluster step,
            // but limit the number of hits.
            if (deltaS < targetLength && curr-first < targetSize) {
                ++curr;
                continue;
            }

            hitIterator hitEnd = *sEnd;
            
            // Check that there are enough hits left for a new cluster.  If
            // not, then quite the loop and add all the remaining points to
            // the last cluster. (not really a good plan, but it's got to
            // suffice for now...
            if (hitEnd-curr < minSize) {
                break;
            }
            
            // Check that the remaining hits cover enough distance along the
            // approximate direction.
            deltaS = dirDiff(*(hitEnd-1), *curr);
            if (deltaS < minLength) break;

            // Check if the not-yet-created cluster could be extended to the
            // end.
            deltaS = dirDiff(*(hitEnd-1), *first);
            if (hitEnd-first < maxSize && deltaS < maxLength) break;
                
            // Time for a new slice of clusters.
            splits.insert(sEnd, curr);

            // Move the the next segment.  This will be the newly created
            // segment.
            sEnd = ++sBegin; ++sEnd;
            
            // Reset first to start looking for a new set of hits.
            first = curr;            
        }

        sEnd = ++sBegin; ++sEnd;
    }
    
    sBegin = splits.begin();
    sEnd = sBegin; ++sEnd;
    while (sEnd != splits.end()) {
        output->push_back(CreateCluster(name,*sBegin,*sEnd,true));
        sBegin = sEnd++;
    }

    return output;
}

template<typename hitIterator>
CP::THandle<CP::TReconObjectContainer> 
CP::CreateShowerClusters(const char* name, hitIterator begin, hitIterator end,
                        const TVector3& approxDir) {
    
    double minLength = 10*unit::mm;
    double maxLength = 30*unit::mm;
    int minSize = std::max(10,(int) (end-begin)/50);
    int maxSize = (end-begin)/5 + 1;
    return CreateClusters(name,begin,end,approxDir,
                          minLength,maxLength,
                          minSize, maxSize);
}


template<typename hitIterator>
CP::THandle<CP::TReconObjectContainer> 
CP::CreateTrackClusters(const char* name, hitIterator begin, hitIterator end,
                        const TVector3& approxDir) {

    double minLength = 5*unit::mm;
    double maxLength = 15*unit::mm;
    int minSize = end - begin;
    minSize = std::max(11, minSize/100);
    int maxSize = end-begin;
    maxSize = maxSize/3 + 1;
    return CreateClusters(name,begin,end,approxDir,
                          minLength,maxLength,
                          minSize, maxSize);
}



#endif
