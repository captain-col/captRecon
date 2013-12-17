#include "TDisassociateHits.hxx"
#include "HitUtilities.hxx"
#include "TPositionDensityCluster.hxx"
#include "CreateCluster.hxx"
#include "CompareReconObjects.hxx"

#include <THandle.hxx>
#include <TReconTrack.hxx>
#include <TReconCluster.hxx>
#include <TCaptLog.hxx>
#include <HEPUnits.hxx>
#include <TUnitsTable.hxx>
#include <TRuntimeParameters.hxx>
 
#include <memory>
#include <set>
#include <cmath>

CP::TDisassociateHits::TDisassociateHits()
    : TAlgorithm("TDisassociateHits", 
                 "Break up objects into separate hits") {
}

CP::TDisassociateHits::~TDisassociateHits() { }

CP::THandle<CP::TAlgorithmResult>
CP::TDisassociateHits::Process(const CP::TAlgorithmResult& input,
                          const CP::TAlgorithmResult&,
                          const CP::TAlgorithmResult&) {
    
    CP::THandle<CP::TReconObjectContainer> inputObjects 
        = input.GetResultsContainer();

    CaptLog("TDisassociateHits Process " << GetEvent().GetContext());

    if (!inputObjects) {
        CaptError("No input objects");
        return CP::THandle<CP::TAlgorithmResult>();
    }

    // Create the output containers.
    CP::THandle<CP::TAlgorithmResult> result = CreateResult();
    std::auto_ptr<CP::TReconObjectContainer> 
        final(new CP::TReconObjectContainer("final"));
    std::auto_ptr<CP::TReconObjectContainer> 
        tracks(new CP::TReconObjectContainer("tracks"));
    std::auto_ptr<CP::TReconObjectContainer> 
        big(new CP::TReconObjectContainer("big"));
    std::auto_ptr<CP::TReconObjectContainer> 
        medium(new CP::TReconObjectContainer("medium"));
    std::auto_ptr<CP::TReconObjectContainer> 
        small(new CP::TReconObjectContainer("small"));

    // Objects to break up.
    std::auto_ptr<CP::TReconObjectContainer> 
        disassociate(new CP::TReconObjectContainer("disassociate"));

    // Divide into objects to save and objects to disintegrate.
    for (CP::TReconObjectContainer::iterator t = inputObjects->begin();
         t != inputObjects->end(); ++t) {
        CP::THandle<CP::TReconTrack> track = *t;
        CP::THandle<CP::TReconCluster> cluster = *t;
        if (track) {
            if (track->GetNodes().size() > 5) {
                // Save all the "big" tracks to final.
                final->push_back(*t);
                tracks->push_back(*t);
                continue;
            }
        }
        else if (!cluster) {
            // Save any other not cluster objects.
            final->push_back(*t);
            continue;
        }
        disassociate->push_back(*t);
    }

    // Get all of the unique hits.
    std::set< CP::THandle<CP::THit> > hits;
    CP::hits::ReconHits(disassociate->begin(), disassociate->end(), hits);

    typedef CP::TPositionDensityCluster< CP::THandle<THit> > ClusterAlgorithm;

    // Find the big clusters.
    ClusterAlgorithm bigClusters(10,20*unit::mm);
    bigClusters.Cluster(hits.begin(), hits.end());

    int nClusters = bigClusters.GetClusterCount();
    CaptNamedLog("TDisassociateHits",
                 "With " << nClusters << " big clusters"
                 << " from " << hits.size() << " hits");
    for (int i=0; i<nClusters; ++i) {
        const ClusterAlgorithm::Points& points 
            = bigClusters.GetCluster(i);
        CP::THandle<CP::TReconCluster> cluster
            = CreateCluster("cluster",points.begin(),points.end());
        big->push_back(cluster);
        final->push_back(cluster);
    }
    
    // Find the medium clusters.
    ClusterAlgorithm midClusters(5,8*unit::mm);
    midClusters.Cluster(bigClusters.GetCluster(nClusters).begin(), 
                        bigClusters.GetCluster(nClusters).end());

    nClusters = midClusters.GetClusterCount();
    CaptNamedLog("TDisassociateHits",
                 "With " << nClusters << " medium clusters"
                 << " from " << hits.size() << " hits");
    for (int i=0; i<nClusters; ++i) {
        const ClusterAlgorithm::Points& points 
            = midClusters.GetCluster(i);
        CP::THandle<CP::TReconCluster> cluster
            = CreateCluster("cluster",points.begin(),points.end());
        medium->push_back(cluster);
        final->push_back(cluster);
    }
    
    // Find the small clusters.
    ClusterAlgorithm smallClusters(1,8*unit::mm);
    smallClusters.Cluster(midClusters.GetCluster(nClusters).begin(), 
                        midClusters.GetCluster(nClusters).end());

    nClusters = smallClusters.GetClusterCount();
    CaptNamedLog("TDisassociateHits",
                 "With " << nClusters << " small clusters"
                 << " from " << hits.size() << " hits");
    for (int i=0; i<nClusters; ++i) {
        const ClusterAlgorithm::Points& points 
            = smallClusters.GetCluster(i);
        CP::THandle<CP::TReconCluster> cluster
            = CreateCluster("cluster",points.begin(),points.end());
        small->push_back(cluster);
        final->push_back(cluster);
    }
    
    std::sort(final->begin(), final->end(), CompareReconObjects());
    std::sort(tracks->begin(), tracks->end(), CompareReconObjects());
    std::sort(big->begin(), big->end(), CompareReconObjects());
    std::sort(medium->begin(), medium->end(), CompareReconObjects());
    std::sort(small->begin(), small->end(), CompareReconObjects());

    result->AddResultsContainer(small.release());
    result->AddResultsContainer(medium.release());
    result->AddResultsContainer(big.release());
    result->AddResultsContainer(tracks.release());
    result->AddResultsContainer(final.release());

    return result;
}
