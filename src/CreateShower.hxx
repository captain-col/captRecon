#ifndef CreateShower_hxx_seen
#define CreateShower_hxx_seen

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
    EXCEPTION(ECreateShower,ECaptRecon);

    /// An exception that there are repeated objects in the track.
    EXCEPTION(EShowerRepeatedObject, ECreateShower);

    /// An exception that there are repeated objects in the track.
    EXCEPTION(EShowerNonCluster, ECreateShower);

    /// Create a vector of clusters that will be used to fill the nodes of the
    /// shower.
    template<typename iterator>
    CP::THandle<CP::TReconObjectContainer> 
    CreateShowerClusters(const char* name, iterator begin, iterator end,
                         const TVector3& approxDir);
    
    /// Create an empty unfitted shower using raw hits.
    template<typename iterator>
    CP::THandle<CP::TReconShower> 
    CreateShower(const char* name, iterator begin, iterator end,
                 const TVector3& approxDir);
    
    /// Create an empty unfit shower using a cluster.
    CP::THandle<CP::TReconShower> 
    CreateShower(const char* name, CP::THandle<CP::TReconCluster>);

    /// Determine the shower width based on a cluster.
    double CreateShowerWidth(CP::THandle<CP::TReconCluster> cluster, 
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
        
    // Set the minimum number hits in a cluster.
    const std::size_t minSize = 5;
    
    // Set the maximum number of hits in a cluster.
    const std::size_t maxSize = std::max((std::size_t) hits.size()/20,
                                         (std::size_t) 10);

    // Set the maximum length of the cluster along the shower direction (after
    // there are minSize hits in the cluster).
    const double minStep = 20*unit::mm;

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

        CP::THandle<CP::TReconCluster> cluster
            = CreateCluster("shower",first,curr);
        output->push_back(cluster);

        // Reset first to start looking for a new set of hits.
        first = curr;
    }

    // Add the last cluster.
    CP::THandle<CP::TReconCluster> cluster = CreateCluster("shower",first,last);
    output->push_back(cluster);
        
    // Find the cluster width and order it so the narrow part is at the front
    // of the container.
    double avg = 0.0;
    double wavg = 0.0;
    double weight = 0.0;
    double index = 0.0;
    for (CP::TReconObjectContainer::iterator c = output->begin();
         c != output->end(); ++c) {
        CP::THandle<CP::TReconCluster> cluster = *c;
        double width = CreateShowerWidth(cluster,dir);
        index += 1;
        avg += index;
        wavg += width*index;
        weight += width;
    }
    avg /= index;
    wavg /= weight;
    if (wavg < avg) std::reverse(output->begin(), output->end());

    return output;
}

template<typename iterator>
CP::THandle<CP::TReconShower> 
CP::CreateShower(const char* name, iterator begin, iterator end,
                 const TVector3& approxDir) {

    CP::THandle<CP::TReconShower> shower(new CP::TReconShower);
    shower->SetAlgorithmName(name);
    shower->SetStatus(CP::TReconBase::kSuccess);
    shower->AddDetector(CP::TReconBase::kTPC);
    shower->SetName("shower");
    
    CP::THitSelection* showerHits = new CP::THitSelection("showerHits");
    std::copy(begin, end, std::back_inserter(*showerHits));
    shower->AddHits(showerHits);
    
    CP::THandle<CP::TReconObjectContainer> clusters
        = CreateShowerClusters(name, 
                               showerHits->begin(), showerHits->end(), 
                               approxDir);

    if (clusters->size() < 5) {
        return CP::THandle<CP::TReconShower>();
    }
    
    // Crudely estimate the direction.
    CP::THandle<CP::TReconCluster> frontCluster = clusters->front();
    CP::THandle<CP::TReconCluster> backCluster = clusters->back();
    
    TVector3 dir = (backCluster->GetPosition().Vect() 
                    - frontCluster->GetPosition().Vect()).Unit();

    TReconNodeContainer& nodes = shower->GetNodes();
    for (CP::TReconObjectContainer::iterator c = clusters->begin();
         c != clusters->end(); ++c) {
        CP::THandle<CP::TReconCluster> cluster = *c;
        CP::THandle<CP::TReconNode> node(new CP::TReconNode);
        CP::THandle<CP::TReconState> state(new CP::TShowerState);
        CP::THandle<CP::TShowerState> showerState = state;
        CP::THandle<CP::TReconBase> object = cluster;
        node->SetState(state);
        node->SetObject(object);
        showerState->SetDirection(dir);
        showerState->SetPosition(cluster->GetPosition());
        showerState->SetEDeposit(cluster->GetEDeposit());
        showerState->SetCone(CreateShowerWidth(cluster,dir));
        nodes.push_back(node);
    }
    
    // Get the state from the front node.  It had better be a showerState!.
    CP::THandle<CP::TShowerState> showerState = shower->GetState();
    CP::THandle<CP::TShowerState> frontState = nodes.front()->GetState();
    *showerState = *frontState;

    return shower;
}

#endif
