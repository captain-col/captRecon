#include "TDisassociateHits.hxx"
#include "HitUtilities.hxx"
#include "TPositionDensityCluster.hxx"
#include "CreateCluster.hxx"
#include "CreateShower.hxx"
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
    std::unique_ptr<CP::TReconObjectContainer> 
        final(new CP::TReconObjectContainer("final"));
    std::unique_ptr<CP::TReconObjectContainer> 
        tracks(new CP::TReconObjectContainer("tracks"));
    std::unique_ptr<CP::TReconObjectContainer>
        showers(new CP::TReconObjectContainer("showers"));
    std::unique_ptr<CP::TReconObjectContainer> 
        big(new CP::TReconObjectContainer("big"));
    std::unique_ptr<CP::TReconObjectContainer> 
        small(new CP::TReconObjectContainer("small"));

    // Check to see if there are tracks and clusters to disassociate.
    std::size_t bigTrackThreshold = 5;
    int trackCount = 0;
    for (CP::TReconObjectContainer::iterator t = inputObjects->begin();
         t != inputObjects->end(); ++t) {
        CP::THandle<CP::TReconTrack> track = *t;
        if (track) ++trackCount;
    }    
    if (trackCount < 15) bigTrackThreshold = 0;

    CaptNamedLog("TDisassociateHits",
                 "Big track threshold " << bigTrackThreshold
                 << " for " << trackCount << " tracks");
    
    // Objects to break up.
    std::unique_ptr<CP::TReconObjectContainer> 
        disassociate(new CP::TReconObjectContainer("disassociate"));

    // Divide into objects to save and objects to disintegrate.
    for (CP::TReconObjectContainer::iterator t = inputObjects->begin();
         t != inputObjects->end(); ++t) {
        CP::THandle<CP::TReconTrack> track = *t;
        CP::THandle<CP::TReconCluster> cluster = *t;
        if (track) {
            if (track->GetNodes().size() > bigTrackThreshold) {
                // Save all the "big" tracks to final.
                final->push_back(*t);
                tracks->push_back(*t);
                continue;
            }
            double length = 0.0;
            for (CP::THitSelection::iterator h = track->GetHits()->begin();
                 h != track->GetHits()->end(); ++h) {
                for (CP::THitSelection::iterator i = h;
                     i != track->GetHits()->end(); ++i) {
                    double v = ((*h)->GetPosition()-(*i)->GetPosition()).Mag();
                    if (v > length) length = v;
                }
            }
            if (length > 100*unit::mm) {
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
    }
    
    // Find the medium clusters.
    ClusterAlgorithm mediumClusters(5,8*unit::mm);
    mediumClusters.Cluster(bigClusters.GetCluster(nClusters).begin(), 
                           bigClusters.GetCluster(nClusters).end());
    
    nClusters = mediumClusters.GetClusterCount();
    CaptNamedLog("TDisassociateHits",
                 "With " << nClusters << " medium clusters"
                 << " from " << bigClusters.GetCluster(nClusters).size() 
                 << " hits");
    
    for (int i=0; i<nClusters; ++i) {
        const ClusterAlgorithm::Points& points 
            = mediumClusters.GetCluster(i);
        CP::THandle<CP::TReconCluster> cluster
            = CreateCluster("cluster",points.begin(),points.end());
        small->push_back(cluster);
    }
    
    // Find the small clusters.
    ClusterAlgorithm smallClusters(1,8*unit::mm);
    smallClusters.Cluster(mediumClusters.GetCluster(nClusters).begin(), 
                            mediumClusters.GetCluster(nClusters).end());
    
    nClusters = smallClusters.GetClusterCount();
    CaptNamedLog("TDisassociateHits",
                 "With " << nClusters << " small clusters"
                 << " from " << mediumClusters.GetCluster(nClusters).size() 
                 << " hits");
    for (int i=0; i<nClusters; ++i) {
        const ClusterAlgorithm::Points& points 
            = smallClusters.GetCluster(i);
        CP::THandle<CP::TReconCluster> cluster
            = CreateCluster("cluster",points.begin(),points.end());
        small->push_back(cluster);
    }
    
    std::sort(final->begin(), final->end(), CompareReconObjects());
    std::sort(tracks->begin(), tracks->end(), CompareReconObjects());
    std::sort(big->begin(), big->end(), CompareReconObjects());
    std::sort(small->begin(), small->end(), CompareReconObjects());

    // Now break up the big clusters as if we were making showers.
    for (CP::TReconObjectContainer::iterator c= big->begin();
         c != big->end(); ++c) {
        CP::THandle<CP::TReconCluster> bigCluster = (*c);
        CP::THandle<CP::TReconShower> shower
            = CreateShower("shower", 
                           bigCluster->GetHits()->begin(),
                           bigCluster->GetHits()->end(),
                           bigCluster->GetLongAxis());
        if (shower) {
            showers->push_back(shower);
            final->push_back(shower);
        }
        else {
            final->push_back(bigCluster);
        }
    }

    result->AddResultsContainer(small.release());
    result->AddResultsContainer(big.release());
    result->AddResultsContainer(showers.release());
    result->AddResultsContainer(tracks.release());
    result->AddResultsContainer(final.release());

    return result;
}
