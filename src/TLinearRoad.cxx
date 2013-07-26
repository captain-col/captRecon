#include "TLinearRoad.hxx"
#include "TClusterMomentsFit.hxx"

#include <HEPUnits.hxx>

#include <TVector3.h>
#include <TPrincipal.h>

#include <numeric>
#include <cmath>
#include <algorithm>
#include <string>
#include <memory>

CP::TLinearRoad::TLinearRoad(const CP::TReconObjectContainer& clusters,
                             const CP::TReconObjectContainer& seed,
                             int maxClusters)
    : fMaxClusters(maxClusters),
      fRoadWidth(12*unit::mm), 
      fRoadStep(10*unit::mm),
      fOpeningAngle(0.15*unit::radian),
      fSeedSize(5),
      fSeedLength(2.0*unit::cm) {

    // Copy the input clusters and make sure the really are all clusters.
    for (CP::TReconObjectContainer::const_iterator c = clusters.begin();
         c != clusters.end(); ++c) {
        CP::THandle<CP::TReconCluster> cluster = *c;
        if (!cluster) {
            CaptError("Non-cluster in input objects");
            std::abort();
        }
        fRemainingClusters.push_back(cluster);
    }
    
    // Copy the seed and make sure it's really all clusters.  This keeps two
    // copies, one (fTrackClusters) will be expanded with new hits, and the
    // other is to keep track of the original seed clusters.
    for (CP::TReconObjectContainer::const_iterator c = seed.begin();
         c != seed.end(); ++c) {
        CP::THandle<CP::TReconCluster> cluster = *c;
        if (!cluster) {
            CaptError("Non-cluster in input objects");
            std::abort();
        }
        fOriginalClusters.push_back(cluster);
        fTrackClusters.push_back(cluster);
    }
}

void CP::TLinearRoad::Process() {
    // Define a queue to hold the current key.
    SeedContainer currentSeed;

    // Find hits at the "upstream" end of the track.
    for (SeedContainer::iterator c = fTrackClusters.begin();
         c != fTrackClusters.end(); ++c) {
        currentSeed.push_back(*c);
        double length = (currentSeed.front()->GetPosition()
                         -currentSeed.back()->GetPosition()).Mag();
        if (currentSeed.size() < fSeedSize) continue;
        if (length < fSeedLength) continue;
        break;
    }

    CaptNamedDebug("road", "Follow road upstream");
    int collectedClusters = 0;
    while (!fRemainingClusters.empty() && currentSeed.size()>2) {
        CP::THandle<CP::TReconCluster> cluster
            = NextCluster(fRemainingClusters, currentSeed, false);

        // If a cluster wasn't found, then stop looking for more.
        if (!cluster) break;

        // Find the cluster being added to the seed and remove it from the
        // remaining clusters.  This keeps the cluster from being found twice.
        RemainingContainer::iterator where =
            std::find(fRemainingClusters.begin(), fRemainingClusters.end(), 
                      cluster);
        fRemainingClusters.erase(where);
        
        // Add the cluster to the seed.  We are searching to the upstream end,
        // so the cluster gets inserted at the beginning.
        currentSeed.push_front(cluster);
        
        while (currentSeed.size() > fSeedSize
               && ((currentSeed.front()->GetPosition()
                    -currentSeed.back()->GetPosition()).Mag()
                   >= fSeedLength)) {
            fTrackClusters.push_front(currentSeed.back());
            currentSeed.pop_back();
        }
        if ((++collectedClusters)>fMaxClusters) break;
    }

    // At the end, flush the clusters in the current seed into the track. 
    std::copy(currentSeed.rbegin(), currentSeed.rend(),
              std::front_inserter(fTrackClusters));
    currentSeed.clear();
    
    // Find clusters at the downstream end of the track.
    for (SeedContainer::reverse_iterator c = fTrackClusters.rbegin();
         c != fTrackClusters.rend();
         ++c) {
        currentSeed.push_front(*c);
        if (currentSeed.size()>=fSeedSize
            && ((currentSeed.front()->GetPosition()
                -currentSeed.back()->GetPosition()).Mag()
                >= fSeedLength)) break;
    }

    collectedClusters = 0;
    CaptNamedDebug("road", "Follow road downstream");
    while (!fRemainingClusters.empty() && currentSeed.size()>2) {
        CP::THandle<CP::TReconCluster> cluster
            = NextCluster(fRemainingClusters, currentSeed, true);

        // If a cluster wasn't found, then stop looking for more.
        if (!cluster) break;

        // Find the cluster being added to the seed and remove it from the
        // remaining clusters.  This keeps the cluster from being found twice.
        RemainingContainer::iterator where 
            = std::find(fRemainingClusters.begin(), fRemainingClusters.end(), 
                        cluster);
        fRemainingClusters.erase(where);
        
        // Add the cluster to the seed.  We are searching to the upstream end,
        // so the cluster gets inserted at the beginning.
        currentSeed.push_back(cluster);
        
        while (currentSeed.size() > fSeedSize
               && ((currentSeed.front()->GetPosition()
                    -currentSeed.back()->GetPosition()).Mag()
                   >= fSeedLength)) {
            fTrackClusters.push_back(currentSeed.front());
            currentSeed.pop_front();
        }
        if ((++collectedClusters)>fMaxClusters) break;
    }
    
    // At the end, flush all the current hits into the track hits 
    std::copy(currentSeed.begin(), currentSeed.end(),
              std::back_inserter(fTrackClusters));

}

CP::THandle<CP::TReconCluster>
CP::TLinearRoad::NextCluster(const RemainingContainer& clusters, 
                             const SeedContainer& seed,
                             bool downstream) {
    
    double bestDistance = 100*unit::meter;
    CP::THandle<CP::TReconCluster> bestCluster;

    TVector3 seedPosition;
    TVector3 seedDirection;

    // Estimate the seed position and direction.
    std::auto_ptr<TPrincipal> principal(new TPrincipal(3,""));
    for (SeedContainer::const_iterator s = seed.begin();s != seed.end(); ++s) {
        double row[3] = {(*s)->GetPosition().X(),
                         (*s)->GetPosition().Y(),
                         (*s)->GetPosition().Z()};
        principal->AddRow(row);
    }
    principal->MakePrincipals();

    // Find the extent of the seed (along the principal axis).
    std::pair<double,double> extent(0,0);
    for (SeedContainer::const_iterator s = seed.begin();s != seed.end(); ++s) {
        for (CP::THitSelection::const_iterator h = (*s)->GetHits()->begin();
             h != (*s)->GetHits()->end(); ++h) {
            double X[3] = {(*h)->GetPosition().X(),
                           (*h)->GetPosition().Y(),
                           (*h)->GetPosition().Z()};
            double P[3] = {0,0,0};
            principal->X2P(X,P);
            if (extent.first > P[0]) extent.first = P[0];
            if (extent.second < P[0]) extent.second = P[0];
        }
    }

    
    
    double minDistance = 100*unit::meter;
    for (RemainingContainer::const_iterator c = clusters.begin();
         c != clusters.end(); ++c) {
        CP::THandle<CP::TReconCluster> cluster = (*c);
        TVector3 diff = cluster->GetPosition().Vect() - seedPosition;
        
        // Only look in the mostly forward direction.  This allows a backward
        // look equal to the road width.
        double dotProd = diff*seedDirection;
        if (dotProd < -0.5*fRoadWidth) {
            CaptNamedDebug("road", "Backward step ");
            continue;
        }

        // The find the local width for a hit at the current distance from the
        // end point.
        double localWidth = fRoadWidth + dotProd*fOpeningAngle;

        // Check that the hit is inside the local width.
        double transDist = (diff - dotProd*seedDirection).Mag();

        // Make sure the hit is inside the road.
        if (transDist > 0.5*localWidth) continue;

        // And keep the hit closest to the current end point.
        if (dotProd < bestDistance) {
            bestCluster = cluster;
            bestDistance = dotProd;
        }
    }
    
    if (bestCluster) {
        CaptNamedDebug("road", "Next cluster at " << minDistance
                        << " at " 
                        << " " << bestCluster->GetPosition().X()
                        << " " << bestCluster->GetPosition().Y()
                        << " " << bestCluster->GetPosition().Z());
    }
    return bestCluster;
}


CP::THandle<CP::TReconTrack> CP::TLinearRoad::GetTrack()  {
    
    CP::THandle<CP::TReconTrack> track(new CP::TReconTrack);
#ifdef BOGUS
    if (fTrackHits.empty()) return track;

    CP::THitSelection* trackHits = new CP::THitSelection("trackHits");
    std::copy(fTrackHits.begin(), fTrackHits.end(),
              std::back_inserter(*trackHits));
    
    track->SetAlgorithmName("TLinearRoad");
    track->SetStatus(CP::TReconBase::kSuccess);
    track->AddDetector(CP::TReconBase::kP0D);
    track->AddHits(trackHits);
    track->SetName(fTrackSeed->GetName());
#endif
    
    return track;     
}
