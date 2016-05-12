#include "TClusterMerge.hxx"
#include "CreateCluster.hxx"
#include "CreateClusters.hxx"
#include "ClusterDistance.hxx"
#include "TPositionNeighbors.hxx"

#include <THandle.hxx>
#include <TReconCluster.hxx>
#include <TCaptLog.hxx>
#include <HEPUnits.hxx>
#include <TUnitsTable.hxx>
#include <TRuntimeParameters.hxx>

#include <TVector3.h>

#include <memory>
#include <cmath>

CP::TClusterMerge::TClusterMerge()
    : TAlgorithm("TClusterMerge", 
                 "Merge clusters that are incorrectly split") {
}

CP::TClusterMerge::~TClusterMerge() { }

CP::THandle<CP::TAlgorithmResult>
CP::TClusterMerge::Process(const CP::TAlgorithmResult& input,
                           const CP::TAlgorithmResult&,
                           const CP::TAlgorithmResult&) {

    CP::THandle<CP::TReconObjectContainer> inputObjects 
        = input.GetResultsContainer();

    CaptLog("TClusterMerge Process " << GetEvent().GetContext());

    if (!inputObjects) {
        CaptError("No input objects");
        return CP::THandle<CP::TAlgorithmResult>();
    }

    // Create the output containers.
    CP::THandle<CP::TAlgorithmResult> result = CreateResult();
    std::unique_ptr<CP::TReconObjectContainer> 
        final(new CP::TReconObjectContainer("final"));

    // Find the clusters in the input object and copy them into the
    // remainingClusters object.  Any non-cluster objects are copied directly
    // to final.
    CP::TReconObjectContainer remainingClusters;
    for (CP::TReconObjectContainer::iterator i = inputObjects->begin();
         i != inputObjects->end(); ++i) {
        CP::THandle<CP::TReconCluster> cluster = *i;
        if (!cluster) {
            final->push_back(*i);
            continue;
        }
        remainingClusters.push_back(*i);
    }

    CP::THandle<CP::TReconCluster> work;
    CP::TReconObjectContainer::iterator next = remainingClusters.begin();
    while (next != remainingClusters.end()) {
        work = *(next++);
        std::cout << "Check next cluster " << remainingClusters.end()-next
                  << std::endl;
        for (CP::TReconObjectContainer::iterator check = next;
             check != remainingClusters.end(); ++check) {
            if (OverlappingClusters(work,*check)) {
                std::cout << "       Combine " << check-next
                          << " " << work->GetHits()->size()
                          << " " << (*check)->GetHits()->size()
                          << std::endl;
                // Move the cluster that will be combined with work off of
                // it's position in the vector.  The cluster to be combined is
                // now at the slot pointed to by next.
                std::swap(*next,*check);
                // Combine the clusters.
                work=CombineClusters(work,*next);
                // After work is combine, start over and check any other
                // clusters to see if they should be combined.  This basically
                // restarts the for loop, but now with one less cluster in the
                // vector to be checked.
                check = ++next;
            }
        }
        TVector3 axis = work->GetLongAxis();
        double width = work->GetMajorAxis().Mag();
        if (width < 0.5*axis.Mag() && axis.Mag() > 15*unit::mm) {
            CP::THandle<CP::TReconObjectContainer> clusters
                = CP::CreateTrackClusters("splitClusters",
                                          work->GetHits()->begin(), 
                                          work->GetHits()->end(),
                                          axis);
            if (clusters) {
                for (CP::TReconObjectContainer::iterator o = clusters->begin();
                     o != clusters->end(); ++o) {
                    final->push_back(*o);
                }
            }
            else {
                final->push_back(work);
            }
        }
        else {
            final->push_back(work);
        }
    };
    
    result->AddResultsContainer(final.release());

    return result;
}

CP::THandle<CP::TReconCluster> CP::TClusterMerge::CombineClusters(
    const CP::THandle<CP::TReconCluster>& cluster1,
    const CP::THandle<CP::TReconCluster>& cluster2) {
    CP::THitSelection hits;
    std::copy(cluster1->GetHits()->begin(),cluster1->GetHits()->end(),
              std::back_inserter(hits));
    std::copy(cluster2->GetHits()->begin(),cluster2->GetHits()->end(),
              std::back_inserter(hits));
    return CreateCluster("Combine",hits.begin(), hits.end());
}

bool CP::TClusterMerge::OverlappingClusters(
    const CP::THandle<CP::TReconCluster>& cluster1,
    const CP::THandle<CP::TReconCluster>& cluster2) {
    double dZ = cluster1->GetPosition().Z() - cluster2->GetPosition().Z();
    const double vicinity = 10*unit::mm;
    // Check if the clusters are close together.
    if (std::abs(dZ) > vicinity) return false;
    if (CP::ClusterVicinity(*cluster1,*cluster2) > vicinity) return false;
    // Check that at least one of the clusters is large.
    const std::size_t minSize = 15; 
    if (cluster1->GetHits()->size()<minSize
        && cluster2->GetHits()->size()<minSize) {
        return false;
    }
    double overlaps = 0;
    for (CP::THitSelection::iterator c1 = cluster1->GetHits()->begin();
         c1 != cluster1->GetHits()->end(); ++c1) {
        for (CP::THitSelection::iterator c2 = cluster2->GetHits()->begin();
             c2 != cluster2->GetHits()->end(); ++c2) {
            double d = std::abs((*c1)->GetPosition().Z()
                                -(*c2)->GetPosition().Z());
            if (d > 5*unit::mm) continue;
            d = std::abs((*c1)->GetPosition().X()-(*c2)->GetPosition().X());
            if (d > 2*unit::mm) continue;
            d = std::abs((*c1)->GetPosition().Y()-(*c2)->GetPosition().Y());
            if (d > 2*unit::mm) continue;
            overlaps += 1.0;
            if (overlaps > 30) return true;
            double r1 = overlaps/cluster1->GetHits()->size();
            double r2 = overlaps/cluster2->GetHits()->size();
            const double minOverlap = 0.1;
            if (r1 > minOverlap) return true;
            if (r2 > minOverlap) return true;
        }
    }
    return false;
}
