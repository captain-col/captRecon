#include "TClusterMerge.hxx"
#include "CreateCluster.hxx"
#include "CreateClusters.hxx"
#include "ClusterDistance.hxx"
#include "TPositionDensityCluster.hxx"

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
    fMinimumClusterSize = 15;
    fClusterSizeFactor = 2.0;
    fClusterOverlap = 0.1;
    fMaximumWidthFraction = 0.25;
    fMinimumLength = 25*unit::mm;
    fClusterExtent = 15*unit::mm;
    fTimeMetric = 0.3;
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
    std::unique_ptr<CP::TReconObjectContainer>
        merged(new CP::TReconObjectContainer("merged"));
    
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
    std::unique_ptr<Neighbors> neighbors;
    while (next != remainingClusters.end()) {
        work = *(next++);
        neighbors.reset(new Neighbors(work->GetHits()->begin(),
                                      work->GetHits()->end()));
        for (CP::TReconObjectContainer::iterator check = next;
             check != remainingClusters.end(); ++check) {
            if (OverlappingClusters(*neighbors,work,*check)) {
                // Move the cluster that will be combined with work off of
                // it's position in the vector.  The cluster to be combined is
                // now at the slot pointed to by next.
                std::swap(*next,*check);
                // Combine the clusters.
                work=CombineClusters(work,*next);
                neighbors.reset(new Neighbors(work->GetHits()->begin(),
                                              work->GetHits()->end()));
                // After work is combine, start over and check any other
                // clusters to see if they should be combined.  This basically
                // restarts the for loop, but now with one less cluster in the
                // vector to be checked.
                check = ++next;
                if (check == remainingClusters.end()) break;
            }
        }

        // Save the result of merging.
        merged->push_back(work);

        // Check if the merged cluster should be split along the long axis.
        // If it should, the split and save the clusters, otherwise, just save
        // the cluster in final.
        TVector3 axis = 2.0*work->GetLongAxis();
        double width = work->GetMajorAxis().Mag();
        if (width < fMaximumWidthFraction*axis.Mag()
            && axis.Mag() > fMinimumLength) {
            // Split the cluster into sub clusters with at least minLength
            // along the long axis.
            double minLength = width;
            // Split with no more than maxLength along the long axis.
            double maxLength = axis.Mag();
            maxLength = std::max(2.0*minLength,maxLength/100.0);
            
            int minSize = work->GetHits()->size();
            minSize = std::max(3, minSize/100);
            int maxSize = work->GetHits()->size();
            maxSize = std::max(4*minSize, maxSize/3);
            
            CP::THandle<CP::TReconObjectContainer> clusters
                = CreateClusters("mergeClusters",
                                 work->GetHits()->begin(), 
                                 work->GetHits()->end(),
                                 axis,
                                 minLength,maxLength,
                                 minSize, maxSize,
                                 false);
            
            if (clusters) {
                for (CP::TReconObjectContainer::iterator o = clusters->begin();
                     o != clusters->end(); ++o) {
                    if ((*o)->GetHits()->size() < 10) {
                        final->push_back(*o);
                        continue;
                    }
                    // Before saving the new split clusters, we need to see if
                    // they should be resplit.  This handles V topologies
                    // where two branches are be merged, but are clearly split
                    // away from the "V".
                    typedef CP::TPositionDensityCluster<CP::THandle<CP::THit> > 
                        ClusterAlgorithm;
                    std::unique_ptr<ClusterAlgorithm> 
                        clusterAlgorithm(
                            new ClusterAlgorithm(1,fClusterExtent));
                    clusterAlgorithm->SetBasis(TVector3(1,0,0),
                                               TVector3(0,1,0),
                                               TVector3(0,0,fTimeMetric));
                    clusterAlgorithm->Cluster((*o)->GetHits()->begin(),
                                              (*o)->GetHits()->end());
                    int nClusters = clusterAlgorithm->GetClusterCount();
                    if (nClusters<2) {
                        final->push_back(*o);
                        continue;
                    }
                    for (int i=0; i<nClusters; ++i) {
                        const ClusterAlgorithm::Points& points 
                            = clusterAlgorithm->GetCluster(i);
                        CP::THandle<CP::TReconCluster> cluster
                            = CreateCluster("mergedCluster",
                                            points.begin(),points.end());
                        final->push_back(cluster);
                    }
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
    
    result->AddResultsContainer(merged.release());
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
    CP::TClusterMerge::Neighbors& neighbors,
    const CP::THandle<CP::TReconCluster>& cluster1,
    const CP::THandle<CP::TReconCluster>& cluster2) {

    // Check if the clusters are close together in Z.
    double dZ = cluster1->GetPosition().Z() - cluster2->GetPosition().Z();
    const double vicinity = 100*unit::mm;
    if (std::abs(dZ) > vicinity) return false;

    // Check that at least one of the clusters is large.
    if (cluster1->GetHits()->size() < (std::size_t) fMinimumClusterSize
        && cluster2->GetHits()->size() < (std::size_t) fMinimumClusterSize) {
        return false;
    }

    double overlaps = 0;
    for (CP::THitSelection::iterator c2 = cluster2->GetHits()->begin();
         c2 != cluster2->GetHits()->end(); ++c2) {

        // Find the closest neighbor in the other cluster.
        Neighbors::iterator neighbor = neighbors.begin((*c2));
        if (neighbor == neighbors.end()) continue;

        // Only look at hits on same triplet of wires.
        double dX = std::abs((*c2)->GetPosition().X()
                     -neighbor->first->GetPosition().X());
        if (dX > 1*unit::mm) continue; 
        double dY = std::abs((*c2)->GetPosition().Y()
                     -neighbor->first->GetPosition().Y());
        if (dY > 1*unit::mm) continue;

        // Find the sizes of the two hits.  If the separating is less than
        // this, they are overlapping.
        double v2 = (*c2)->GetUncertainty().Z();
        double v1 = neighbor->first->GetUncertainty().Z();
        double v = fClusterSizeFactor*(std::abs(v1)+std::abs(v2));

        // The hit size is always bigger than 1.5*mm (total of 3 mm for two
        // hits) since that the resolution scale of our detectors.
        v = std::max(v,3.0*unit::mm);

        // Find the distance between the hits in the Z direction.
        double dZ = std::abs((*c2)->GetPosition().Z()
                            -neighbor->first->GetPosition().Z());
        if (dZ > v) continue;
        
        // We have an overlap.
        overlaps += 1.0;

        // If we have more than 30 overlaps, then the two clusters should be
        // combined.  This speeds things up.
        if (overlaps > 30) return true;

        double r1 = overlaps/cluster1->GetHits()->size();
        double r2 = overlaps/cluster2->GetHits()->size();
        if (r1 > fClusterOverlap) return true;
        if (r2 > fClusterOverlap) return true;
    }
    
    return false;
}
