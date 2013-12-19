#ifndef CreateShower_hxx_seen
#define CreateShower_hxx_seen

#include <TReconCluster.hxx>
#include <TReconShower.hxx>
#include <THandle.hxx>
#include <THit.hxx>
#include <HEPUnits.hxx>

#include <algorithm>

namespace CP {

    /// A base exception for the create track template.
    EXCEPTION(ECreateTrack,ECaptRecon);

    /// An exception that there are repeated objects in the track.
    EXCEPTION(ETrackRepeatedObject, ECreateTrack);

    /// An exception that there are repeated objects in the track.
    EXCEPTION(ETrackNonCluster, ECreateTrack);

    template<typename iterator>
    CP::THandle<CP::TReconObjectContainer> 
    CreateShowerClusters(const char* name, iterator begin, iterator end,
                         const TVector3& approxDir);

};

namespace {
    struct CreateShowerDirCompare {
        explicit CreateShowerDirCompare(const TVector3& dir) 
            : fDir(dir) {}
        bool operator() (const CP::THandle<CP::THit> lhs,
                         const CP::THandle<CP::THit> rhs) {
            return fDir*lhs->GetPosition() < fDir*rhs->GetPosition();
            
        }
    private:
        TVector3 fDir;
    };
};
    
/// Take iterators to a bunch of hits and slice them up into clusters along
/// the approximate direction.
template<typename iterator>
CP::THandle<CP::TReconObjectContainer> 
CP::CreateShowerClusters(const char* name, iterator begin, iterator end,
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
    std::sort(hits.begin(), hits.end(), CreateShowerDirCompare(dir));
        
    // Create clusters in order of the hits.
    std::size_t minSize = 10;
    std::size_t maxSize = 0.02*hits.size();
    maxSize = std::max((std::size_t) 15,maxSize);

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
            
        // Make sure the cluster covers a minimum range in Z (determined by
        // hit resolution).
        double minStep = 20*unit::mm;

        // Make sure that the cluster at least the expected cluster step, but
        // limit the number of hits.
        if (deltaS < minStep && curr-first < maxSize) {
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

        std::cout << "Create " << curr-first 
                  << " " << deltaS
                  << " " << minStep 
                  << " " << maxSize
                  << std::endl;

        CP::THandle<CP::TReconCluster> cluster
            = CreateCluster("shower",first,curr);
        output->push_back(cluster);

        // Reset first to start looking for a new set of hits.
        first = curr;
    }

    // Add the last cluster.
    CP::THandle<CP::TReconCluster> cluster = CreateCluster("shower",first,last);
    output->push_back(cluster);
        
    return output;
}

#endif
