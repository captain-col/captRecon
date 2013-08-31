#include "TLinearRoad.hxx"
#include "TClusterMomentsFit.hxx"
#include "HitUtilities.hxx"

#include <ostreamTLorentzVector.hxx>
#include <ostreamTVector3.hxx>

#include <HEPUnits.hxx>

#include <TVector3.h>
#include <TPrincipal.h>

#include <numeric>
#include <cmath>
#include <algorithm>
#include <string>
#include <memory>

bool 
CheckUniqueInternal(std::vector< CP::THandle<CP::TReconCluster> >& vec) {
    std::sort(vec.begin(),vec.end());
    std::vector< CP::THandle<CP::TReconCluster> >::iterator 
        end = std::unique(vec.begin(), vec.end());
    if (end != vec.end()) {
        CaptLog("Duplicate objects in container");
        return false;
    }
    return true;
}

template <typename iterator>
bool CheckUnique(iterator begin, iterator end) {
    std::vector< CP::THandle<CP::TReconCluster> > vec;
    std::copy(begin, end, std::back_inserter(vec));
    if (!CheckUniqueInternal(vec)) {
        for (iterator i = begin; i != end; ++i) {
            CaptLog("   " << CP::GetPointer(*i));
        }
        return false;
    }
    return true;
}

CP::TLinearRoad::TLinearRoad(int maxClusters)
    : fMaxClusters(maxClusters),
      fRoadWidth(12*unit::mm), 
      fRoadStep(10*unit::mm),
      fOpeningAngle(0.15*unit::radian),
      fSeedSize(5),
      fSeedLength(2.0*unit::cm) { }

void CP::TLinearRoad::FillRemains(CP::TReconObjectContainer& remains) const {
    remains.clear();
    std::copy(fRemainingClusters.begin(), fRemainingClusters.end(),
              std::back_inserter(remains));
}

void CP::TLinearRoad::Process(const CP::TReconObjectContainer& seed,
                              const CP::TReconObjectContainer& clusters) {

    CaptLog("Follow road from " << seed.size() << " cluster seed"
            << " with " << clusters.size() << " remaining clusters.");

    // Check the inputs and internal structures to make sure things are OK.
    if (!fRemainingClusters.empty()) {
        CaptError("Remaining clusters not empty.");
        std::abort();
    }

    if (!fOriginalClusters.empty()) {
        CaptError("Original clusters not empty.");
        std::abort();
    }

    if (!fTrackClusters.empty()) {
        CaptError("Track clusters not empty.");
        std::abort();
    }

    // Copy the input clusters and make sure they really are all clusters.
    fRemainingClusters.clear();
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

    // Define a queue to hold the current key.
    SeedContainer currentSeed;

    // Find hits at the "upstream" end of the track.
    for (SeedContainer::iterator c = fTrackClusters.begin();
         c != fTrackClusters.end(); ++c) {
        currentSeed.push_back(*c);
        double length = (currentSeed.front()->GetPosition().Vect()
                         -currentSeed.back()->GetPosition().Vect()).Mag();
        CaptNamedDebug("road", "Add upstream " << 
                       currentSeed.back()->GetPosition().Vect()
                       << "   size: " << currentSeed.size() 
                       << "   length: " << length);
        if (currentSeed.size() < fSeedSize) continue;
        if (length < fSeedLength) continue;
        break;
    }

    for (SeedContainer::iterator s = currentSeed.begin(); 
         s != currentSeed.end(); ++s) {
        SeedContainer::iterator where =
            std::find(fTrackClusters.begin(), fTrackClusters.end(), 
                      *s);
        if (where ==  fTrackClusters.end()) {
            CaptError("Seed not found in track.");
            continue;
        }
        fTrackClusters.erase(where);
    }

    CaptNamedDebug("road", "Follow road upstream");
    int collectedClusters = 0;
    while (!fRemainingClusters.empty() && currentSeed.size()>2) {
        CP::THandle<CP::TReconCluster> cluster
            = NextCluster(currentSeed, fRemainingClusters, false);

        // If a cluster wasn't found, then stop looking for more.
        if (!cluster) break;
    
        // Find the cluster being added to the seed and remove it from the
        // remaining clusters.  This keeps the cluster from being found twice.
        RemainingContainer::iterator where =
            std::find(fRemainingClusters.begin(), fRemainingClusters.end(), 
                      cluster);
        if (where ==  fRemainingClusters.end()) {
            CaptError("Upstream cluster not found in remaining clusters");
        }
        fRemainingClusters.erase(where);
        
        // Add the cluster to the seed.  We are searching to the upstream end,
        // so the cluster gets inserted at the beginning.
        currentSeed.push_front(cluster);
        
        while (currentSeed.size() > fSeedSize
               && ((currentSeed.front()->GetPosition().Vect()
                    -currentSeed.back()->GetPosition().Vect()).Mag()
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
        double length = (currentSeed.front()->GetPosition().Vect()
                         -currentSeed.back()->GetPosition().Vect()).Mag();
        CaptNamedDebug("road", "Add downstream " << 
                       currentSeed.front()->GetPosition().Vect()
                       << "   size: " << currentSeed.size() 
                       << "   length: " << length);
        if (currentSeed.size() < fSeedSize) continue;
        if (length < fSeedLength) continue;
        break;
    }

    for (SeedContainer::iterator s = currentSeed.begin(); 
         s != currentSeed.end(); ++s) {
        SeedContainer::iterator where =
            std::find(fTrackClusters.begin(), fTrackClusters.end(), 
                      *s);
        if (where ==  fTrackClusters.end()) {
            CaptError("Seed not found in track.");
            continue;
        }
        fTrackClusters.erase(where);
    }

    collectedClusters = 0;
    CaptNamedDebug("road", "Follow road downstream");
    while (!fRemainingClusters.empty() && currentSeed.size()>2) {
        CP::THandle<CP::TReconCluster> cluster
            = NextCluster(currentSeed, fRemainingClusters, true);

        // If a cluster wasn't found, then stop looking for more.
        if (!cluster) break;

        // Find the cluster being added to the seed and remove it from the
        // remaining clusters.  This keeps the cluster from being found twice.
        RemainingContainer::iterator where 
            = std::find(fRemainingClusters.begin(), fRemainingClusters.end(), 
                        cluster);
        if (where == fRemainingClusters.end()) {
            CaptError("Downstream cluster not found in remaining clusters");
        }
        fRemainingClusters.erase(where);
        
        // Add the cluster to the seed.  We are searching to the upstream end,
        // so the cluster gets inserted at the beginning.
        currentSeed.push_back(cluster);
        
        while (currentSeed.size() > fSeedSize
               && ((currentSeed.front()->GetPosition().Vect()
                    -currentSeed.back()->GetPosition().Vect()).Mag()
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

double 
CP::TLinearRoad::FindPositionPrincipal(TPrincipal& pca,
                                       const TVector3& position) {
    double X[3] = {position.X(), position.Y(), position.Z()};
    double P[3] = {0,0,0};
    pca.X2P(X,P);
    return P[0];
}

TVector3 
CP::TLinearRoad::FindPrincipalPosition(TPrincipal& pca,
                                       double principal) {
    double X[3] = {0,0,0};
    double P[3] = {principal,0,0};
    pca.P2X(P,X,3);
    return TVector3(X[0],X[1],X[2]);
}

CP::THandle<CP::TReconCluster>
CP::TLinearRoad::NextCluster(const SeedContainer& seed,
                             const RemainingContainer& clusters, 
                             bool extendBack) {
    
    double bestDistance = 100*unit::meter;
    CP::THandle<CP::TReconCluster> bestCluster;

    // The direction is estimated using a PCA analysis of the seed cluster
    // positions.
    std::auto_ptr<TPrincipal> principal(new TPrincipal(3,""));
    for (SeedContainer::const_iterator s = seed.begin();s != seed.end(); ++s) {
        double row[3] = {(*s)->GetPosition().X(),
                         (*s)->GetPosition().Y(),
                         (*s)->GetPosition().Z()};
        principal->AddRow(row);
    }
    principal->MakePrincipals();

    // Find the extent of the seed (along the principal axis).  The extent is
    // calculated based on the min and max hit positions along the major axis.
    // The extent is then used to find the seedPosition.
    std::pair<double,double> extent(0,0);
    for (SeedContainer::const_iterator s = seed.begin();s != seed.end(); ++s) {
        for (CP::THitSelection::const_iterator h = (*s)->GetHits()->begin();
             h != (*s)->GetHits()->end(); ++h) {
            double p = FindPositionPrincipal(*principal,
                                             (*h)->GetPosition());
            if (extent.first > p) extent.first = p;
            if (extent.second < p) extent.second = p;
        }
    }

    // Find the position along the principal axis for the front and back of
    // the seed and use that to determine the end of the seed to be extended.
    double pFront = FindPositionPrincipal(*principal,
                                          seed.front()->GetPosition().Vect());
    double pBack = FindPositionPrincipal(*principal,
                                         seed.back()->GetPosition().Vect());
    double pPosition;

    if (extendBack) {
        if (pFront < pBack) pPosition = extent.second;
        else pPosition = extent.first;
    }
    else {
        if (pFront < pBack) pPosition = extent.first;
        else pPosition = extent.second;
    }
    TVector3 seedPosition = FindPrincipalPosition(*principal,pPosition);
    TVector3 seedDirection = FindPrincipalPosition(*principal,0.0);
    seedDirection = (seedPosition - seedDirection).Unit();

    // Check all of the clusters to see which one gets added next.
    for (RemainingContainer::const_iterator c = clusters.begin();
         c != clusters.end(); ++c) {
        CP::THandle<CP::TReconCluster> cluster = (*c);
        TVector3 diff = cluster->GetPosition().Vect() - seedPosition;
        
        // Only look in the mostly forward direction.  This allows a backward
        // look equal to the road width.
        double dotProd = diff*seedDirection;
        if (dotProd < -0.5*fRoadWidth) {
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
        CaptNamedDebug("road", "Next cluster at " << bestDistance
                       << " at " << bestCluster->GetPosition().Vect());
    }
    return bestCluster;
}

// This method creates a TTrackState from a position and direction (with
// covariance matrices) and a THit.
CP::THandle<CP::TTrackState> 
CP::TLinearRoad::CreateTrackState(CP::THandle<CP::TReconCluster> object,
                                  const TVector3& direction) {
    CP::THandle<CP::TTrackState> tstate(new CP::TTrackState);
    
    // Set the EDeposit
    tstate->SetEDeposit(object->GetEDeposit());
    tstate->SetEDepositVariance(object->GetEDeposit());
    
    // Set the value and covariance matrix for the position
    tstate->SetPosition(object->GetPosition());      
    tstate->SetPositionVariance(object->GetPositionVariance().X(),
                                object->GetPositionVariance().Y(),
                                object->GetPositionVariance().Z(),
                                object->GetPositionVariance().T());
    
    tstate->SetDirection(direction);
    tstate->SetDirectionVariance(0,0,0);
    
    // Set the curvature.  It should be zero and fixed, since we assume
    // straight line.
    tstate->SetCurvature(0.0);
    tstate->SetCurvatureVariance(0.0);
    
    return tstate;
}

CP::THandle<CP::TReconTrack> CP::TLinearRoad::GetTrack()  {
    CP::THandle<CP::TReconTrack> track(new CP::TReconTrack);
    if (fTrackClusters.size() < 2) return track;

    CheckUnique(fTrackClusters.begin(), fTrackClusters.end());

    track->SetAlgorithmName("TLinearRoad");
    track->SetStatus(CP::TReconBase::kSuccess);
    track->AddDetector(CP::TReconBase::kTPC);
    track->SetName("track");

    CP::THandle<CP::THitSelection> hits 
        = CP::hits::ReconHits(fTrackClusters.begin(),fTrackClusters.end());
    CP::THitSelection* trackHits = new CP::THitSelection("trackHits");
    std::copy(hits->begin(), hits->end(),std::back_inserter(*trackHits));
    track->AddHits(trackHits);
    
    TReconNodeContainer& nodes = track->GetNodes();
    for (SeedContainer::iterator c = fTrackClusters.begin();
         c != fTrackClusters.end(); ++c) {
        CP::THandle<CP::TReconNode> node(new CP::TReconNode);
        
        TVector3 dir;
        if ((c+1) != fTrackClusters.end()) {
            dir = (*(c+1))->GetPosition().Vect()-(*c)->GetPosition().Vect();
        }
        else {
            dir = (*c)->GetPosition().Vect()-(*(c-1))->GetPosition().Vect();
        }
        dir = dir.Unit();
	CP::THandle<CP::TTrackState> nodeState = CreateTrackState((*c), dir);
        
        CP::THandle<CP::TReconState> castState = nodeState;
        node->SetState(castState);
        CP::THandle<CP::TReconBase> castObject = *c;
        node->SetObject(castObject);
        nodes.push_back(node);
    }


    return track;     
}
