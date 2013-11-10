#include "TClusterSlice.hxx"
#include "HitUtilities.hxx"
#include "TLinearRoad.hxx"
#include "TBestTube.hxx"
#include "TTrackFit.hxx"
#include "THitProximityCluster.hxx"

#include <THandle.hxx>
#include <TReconCluster.hxx>
#include <TCaptLog.hxx>
#include <HEPUnits.hxx>
#include <TRuntimeParameters.hxx>

#include <memory>
#include <iostream>
#include <sstream>
#include <cmath>

namespace {
    struct HitZSort {
        bool operator() (CP::THandle<CP::THit> lhs, 
                           CP::THandle<CP::THit> rhs) {
            return (lhs->GetPosition().Z() < rhs->GetPosition().Z());
        }
    };
};


CP::TClusterSlice::TClusterSlice()
    : TAlgorithm("TClusterSlice", "Find Simply Connected Hits") {
    fMinPoints = CP::TRuntimeParameters::Get().GetParameterI(
        "captRecon.clusterSlice.minPoints");
    fMaxDist = CP::TRuntimeParameters::Get().GetParameterD(
        "captRecon.clusterSlice.maxDistance");
    fMinHits = CP::TRuntimeParameters::Get().GetParameterD(
        "captRecon.clusterSlice.minHits");
    fClusterStep = CP::TRuntimeParameters::Get().GetParameterD(
        "captRecon.clusterSlice.clusterStep");
}

CP::TClusterSlice::~TClusterSlice() { }

CP::THandle<CP::TReconObjectContainer> 
CP::TClusterSlice::MakeSlices(CP::THandle<CP::THitSelection> inputHits) {
    CP::THandle<CP::TReconObjectContainer> result(
        new CP::TReconObjectContainer("result"));

    if (!inputHits) return result;
    if (inputHits->size() < fMinHits) return result;
    
    /// Copy the input hits and sort them by the Z position.
    std::auto_ptr<CP::THitSelection> hits(new THitSelection);
    hits->reserve(inputHits->size());
    std::copy(inputHits->begin(), inputHits->end(), std::back_inserter(*hits));
    std::sort(hits->begin(), hits->end(), HitZSort());
    
    // Distance between first and last point.
    double deltaZ = std::abs(hits->front()->GetPosition().Z() 
                             - hits->back()->GetPosition().Z());
    // The number of steps based the cluster step in Z.
    int steps = 1 + deltaZ/fClusterStep;
    // The adjusted step size so that each step is the same "thickness" in Z.
    double zStep = deltaZ / steps;

    CaptNamedLog("TClusterSlice","Make " << steps << " slices"
                 << " in " << deltaZ << " mm"
                 << " (" << zStep << " mm each)");

    std::auto_ptr<CP::HitProximity::Cluster> 
        clusterAlgorithm(new CP::HitProximity::Cluster(1,6*unit::mm));

    CP::THitSelection::iterator curr = hits->begin();
    CP::THitSelection::iterator end = hits->end();
    CP::THitSelection::iterator first = curr;
    while (curr != end) {
        // Make sure each cluster has at least 2*fMinPoints)
        if (curr-first < 2*fMinPoints) {
            ++curr;
            continue;
        }
        // Make sure that the cluster covers a range in Z.
        deltaZ = std::abs((*curr)->GetPosition().Z()
                          -(*first)->GetPosition().Z());
        if (deltaZ < zStep) {
            ++curr;
            continue;
        }
        // Time for a new cluster.
        ++curr;
        // Check that there are enough hits left for a new cluster.  If not,
        // then add all the remaining points to this cluster. (not really a
        // good plan, but it's got to suffice for now...}
        if (end-curr < 2*fMinPoints) break;
        // The hits between first and curr should be run through the density
        // cluster again since it's very likely that the hits are disjoint in
        // the Z slice.
        clusterAlgorithm->Cluster(first,curr);
        int nClusters = clusterAlgorithm->GetClusterCount();
        CaptNamedVerbose("TClusterSlice","    Slice with " << nClusters
                     << " clusters from " << curr-first << " hits");
        for (int i=0; i<nClusters; ++i) {
            const CP::HitProximity::Cluster::Points& points 
                = clusterAlgorithm->GetCluster(i);
            CaptNamedVerbose("TClusterSlice","       Cluster " << i
                         << " with " << points.size() << " hits");
            CP::THandle<CP::TReconCluster> cluster(new CP::TReconCluster);
            cluster->FillFromHits("zCluster",points.begin(),points.end());
            result->push_back(cluster);
        }
        // Reset first to start looking for a new set of hits.
        first = curr;
    }
    if (first != end) {
        // Build the final clusters.
        clusterAlgorithm->Cluster(first,end);
        int nClusters = clusterAlgorithm->GetClusterCount();
        CaptNamedVerbose("TClusterSlice","    Final slice with " << nClusters
                     << " clusters from " << end-first << " hits");
        for (int i=0; i<nClusters; ++i) {
            const CP::HitProximity::Cluster::Points& points 
                = clusterAlgorithm->GetCluster(i);
            CaptNamedVerbose("TClusterSlice","       Cluster " << i
                         << " with " << points.size() << " hits");
            CP::THandle<CP::TReconCluster> cluster(new CP::TReconCluster);
            cluster->FillFromHits("zCluster",points.begin(),points.end());
            result->push_back(cluster);
        }
    }

    return result;
}

CP::THandle<CP::TReconObjectContainer>
CP::TClusterSlice::MakeTracks(const CP::TReconObjectContainer& input,
                              CP::THandle<CP::TAlgorithmResult> thisResult) {
    CP::THandle<CP::TReconObjectContainer> 
        output(new CP::TReconObjectContainer);
    CP::TReconObjectContainer remains;
    CP::TReconObjectContainer seed;
    std::copy(input.begin(), input.end(), std::back_inserter(remains));
    CP::TTrackFit trackFitter;
    for (int i=0; i<50; ++i) {
        CaptNamedLog("TClusterSlice","Look for seed " << i
                     << " in " << remains.size() << " clusters");
        // Find a seed.
        seed.clear();
        std::auto_ptr<CP::TBestTube> bestTube(new CP::TBestTube);
        std::size_t lastRemains = remains.size();
        bestTube->Process(remains);
        bestTube->FillSeed(seed);
        bestTube->FillRemains(remains);
        CaptNamedLog("TClusterSlice","Found seed with " 
                     << seed.size() << " clusters");
        if (seed.size()<3) break;

#ifdef DEBUG_SEED_CLUSTERS
        // If this is for debugging, save the seed clusters.
        if (thisResult) {
            std::ostringstream seedName;
            seedName << "seed" << i;
            std::auto_ptr<CP::TReconObjectContainer> 
                saveSeed(new CP::TReconObjectContainer(seedName.str().c_str()));
            std::copy(seed.begin(), seed.end(), std::back_inserter(*saveSeed));
            thisResult->AddResultsContainer(saveSeed.release());
        }
#endif

        // Make the track.
        std::auto_ptr<CP::TLinearRoad> road(new CP::TLinearRoad);
        road->Process(seed,remains);
        CP::THandle<CP::TReconTrack> track = road->GetTrack();

#ifdef DEBUG_TRACK_CLUSTERS
        // If this is for debugging, save the track clusters.
        if (thisResult) {
            std::ostringstream trackName;
            trackName << "rawTrack" << i;
            std::auto_ptr<CP::TReconObjectContainer> 
                saveTrack(new CP::TReconObjectContainer(
                              trackName.str().c_str()));
            for (CP::TReconNodeContainer::iterator n 
                     = track->GetNodes().begin();
                 n != track->GetNodes().end(); ++n) {
                saveTrack->push_back((*n)->GetObject());
            }
            thisResult->AddResultsContainer(saveTrack.release());
        }
#endif

        // Get the remaining hits.
        road->FillRemains(remains);
        if (track) {
            CaptNamedLog("TClusterSlice","Track found with "
                         << track->GetNodes().size() << " nodes"
                         << " and " << remains.size() << " remaining clusters");
            // Fit the track and then save it to the output.
            track = trackFitter(track);
            output->push_back(track);
        }
        // Check if we continue.
        if (remains.size()<2 || lastRemains == remains.size()) break;
    }
    std::copy(remains.begin(), remains.end(), std::back_inserter(*output));
    return output;
}

CP::THandle<CP::TAlgorithmResult>
CP::TClusterSlice::Process(const CP::TAlgorithmResult& input,
                           const CP::TAlgorithmResult&,
                           const CP::TAlgorithmResult&) {

    CP::THandle<CP::TReconObjectContainer> inputObjects 
        = input.GetResultsContainer();
    CP::THandle<CP::THitSelection> inputHits = input.GetHits();

    if (!inputObjects && !inputHits) {
        CaptError("No input objects or hits");
        return CP::THandle<CP::TAlgorithmResult>();
    }

    CaptLog("TClusterSlice Process " << GetEvent().GetContext());
    
    // Create the output containers.
    CP::THandle<CP::TAlgorithmResult> result = CreateResult();
    std::auto_ptr<CP::TReconObjectContainer> 
        final(new CP::TReconObjectContainer("final"));
    std::auto_ptr<CP::TReconObjectContainer> 
        slices(new CP::TReconObjectContainer("slices"));
    std::auto_ptr<CP::THitSelection> used(new CP::THitSelection("used"));

    if (inputObjects) {
        CaptLog("   Using " <<  inputObjects->size() << " objects");
        for (CP::TReconObjectContainer::iterator obj = inputObjects->begin();
             obj != inputObjects->end(); ++obj) {
            CP::THandle<CP::THitSelection> hits = (*obj)->GetHits(); 
            if (!hits) {
                CaptError("Input object without hits");
                continue;
            }
            CaptLog("   Hits in object: " << hits->size());
            CP::THandle<CP::TReconObjectContainer> cuts = MakeSlices(hits);
            CaptLog("   Slices made: " << cuts->size());
            std::copy(cuts->begin(), cuts->end(), std::back_inserter(*slices));
            CP::THandle<CP::TReconObjectContainer> tracks 
                = MakeTracks(*cuts,result);
            CaptLog("   Tracks found: " << tracks->size());
            std::copy(tracks->begin(),tracks->end(),std::back_inserter(*final));
        }
    }
    else if (inputHits) {
        CaptLog("   Using " <<  inputHits->size() << " hits");
        CP::THandle<CP::TReconObjectContainer> cuts = MakeSlices(inputHits);
        CaptLog("   Slices made: " << cuts->size());
        std::copy(cuts->begin(), cuts->end(), std::back_inserter(*slices));
        CP::THandle<CP::TReconObjectContainer> tracks = MakeTracks(*cuts,
            result);
        CaptLog("   Tracks found: " << tracks->size());
        std::copy(tracks->begin(),tracks->end(),std::back_inserter(*final));
    }


    // Copy all of the hits that got added to a reconstruction object into the
    // used hit selection.
    CP::THandle<CP::THitSelection> hits 
        = CP::hits::ReconHits(final->begin(), final->end());
    if (hits) {
        used->reserve(hits->size());
        std::copy(hits->begin(), hits->end(), std::back_inserter(*used));
    }

    result->AddHits(used.release());
    result->AddResultsContainer(slices.release());
    result->AddResultsContainer(final.release());

    return result;
}
