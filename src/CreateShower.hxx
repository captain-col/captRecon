#ifndef CreateShower_hxx_seen
#define CreateShower_hxx_seen

#include "CreateClusters.hxx"

#include <ECore.hxx>
#include <TReconCluster.hxx>
#include <TReconShower.hxx>
#include <THandle.hxx>
#include <THit.hxx>

#include <algorithm>

namespace CP {

    /// A base exception for the create track template.
    EXCEPTION(ECreateShower,ECaptRecon);

    /// An exception that there are repeated objects in the track.
    EXCEPTION(EShowerRepeatedObject, ECreateShower);

    /// An exception that there are repeated objects in the track.
    EXCEPTION(EShowerNonCluster, ECreateShower);

    /// Create an empty unfitted shower using raw hits.
    template<typename hitIterator>
    CP::THandle<CP::TReconShower> 
    CreateShower(const char* name, hitIterator begin, hitIterator end,
                 const TVector3& approxDir);
    
    /// Create an empty unfit shower using a cluster.
    CP::THandle<CP::TReconShower> 
    CreateShower(const char* name, CP::THandle<CP::TReconCluster>);

    /// Determine the shower width based on a cluster.
    double CreateShowerWidth(CP::THandle<CP::TReconCluster> cluster, 
                             const TVector3& approxDir);
};

//////////////////////////////////////////////////////////////////
// IMPLEMENTATION
//////////////////////////////////////////////////////////////////

template<typename hitIterator>
CP::THandle<CP::TReconShower> 
CP::CreateShower(const char* name, hitIterator begin, hitIterator end,
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
        = CreateShowerClusters(name, showerHits->begin(), showerHits->end(), 
                               approxDir);
    
    // Find the cluster width and order it so the narrow part is at the front
    // of the container.
    double avg = 0.0;
    double wavg = 0.0;
    double weight = 0.0;
    double index = 0.0;
    for (CP::TReconObjectContainer::iterator c = clusters->begin();
         c != clusters->end(); ++c) {
        CP::THandle<CP::TReconCluster> cluster = *c;
        double width = CreateShowerWidth(cluster,approxDir);
        index += 1;
        avg += index;
        wavg += width*index;
        weight += width;
    }
    avg /= index;
    wavg /= weight;
    if (wavg < avg) std::reverse(clusters->begin(), clusters->end());

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
