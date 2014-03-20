#ifndef CreateClusters_hxx_seen
#define CreateClusters_hxx_seen

#include "CreateCluster.hxx"

#include <ECore.hxx>
#include <TReconCluster.hxx>
#include <TReconShower.hxx>
#include <THandle.hxx>
#include <THit.hxx>
#include <HEPUnits.hxx>

#include <algorithm>

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
                   const TVector3& approxDir);
    
};

//////////////////////////////////////////////////////////////////
// IMPLEMENTATION
//////////////////////////////////////////////////////////////////

namespace {
    struct CreateClustersDirCompare {
        explicit CreateClustersDirCompare(const TVector3& dir) 
            : fDir(dir) {}
        bool operator() (const CP::THandle<CP::THit> lhs,
                         const CP::THandle<CP::THit> rhs) {
            return fDir*lhs->GetPosition() < fDir*rhs->GetPosition();
            
        }
    private:
        TVector3 fDir;
    };
};
    
template<typename hitIterator>
CP::THandle<CP::TReconObjectContainer> 
CP::CreateClusters(const char* name, hitIterator begin, hitIterator end,
                         const TVector3& approxDir) {

    TVector3 dir = approxDir;
    dir = dir.Unit();

    // Make the output object.
    CP::THandle<CP::TReconObjectContainer> 
        output(new CP::TReconObjectContainer(name));

    // copy the hits into local storage so they can be sorted.
    CP::THitSelection hits;
    std::copy(begin, end,std::back_inserter(hits));

    // Sort the hits along the direction of the shower.
    std::sort(hits.begin(), hits.end(), CreateClustersDirCompare(dir));
        
    // Set the minimum number hits in a cluster.  A cluster must have at least
    // this many hits.
    const std::size_t minSize = std::max((std::size_t) hits.size()/50,
                                         (std::size_t) 10);
    
    // Set the maximum number of hits in a cluster.
    const std::size_t maxSize = std::max((std::size_t) hits.size()/20,
                                         (std::size_t) 10);

    // Set the maximum length of the cluster along the shower direction (after
    // there are minSize hits in the cluster).
    const double maxStep = 20*unit::mm;

    // Create clusters in order of the hits.
    CP::THitSelection::iterator curr = hits.begin();
    CP::THitSelection::iterator last = hits.end();
    CP::THitSelection::iterator first = curr;
    while (curr != last) {

        // Make sure each cluster has at least clustSize hits.
        if (curr-first < minSize) {
            ++curr;
            continue;
        }
            
        // Make sure the cluster isn't too long.
        double deltaS = dir*(*curr)->GetPosition()
            -dir*(*first)->GetPosition();
            
        // Make sure that the cluster at least the expected cluster step, but
        // limit the number of hits.
        if (deltaS < maxStep && curr-first < maxSize) {
            ++curr;
            continue;
        }
            
        // Time for a new slice of clusters.
        ++curr;
            
        // Check that there are enough hits left for a new cluster.  If not,
        // then quite the loop and add all the remaining points to the last
        // cluster. (not really a good plan, but it's got to suffice for
        // now...
        if (last-curr < minSize) break;

        CP::THandle<CP::TReconCluster> cluster
            = CreateCluster("shower",first,curr,true);
        output->push_back(cluster);

        // Reset first to start looking for a new set of hits.
        first = curr;
    }

    // Add the last cluster.
    CP::THandle<CP::TReconCluster> cluster 
        = CreateCluster("clusters",first,last);
    output->push_back(cluster);
        
    return output;
}
#endif
