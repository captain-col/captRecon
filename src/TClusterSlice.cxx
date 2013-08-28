#include "TClusterSlice.hxx"
#include "TTmplDensityCluster.hxx"
#include "HitUtilities.hxx"
#include "TLinearRoad.hxx"
#include "TBestTube.hxx"
#include "TTrackFit.hxx"

#include <THandle.hxx>
#include <TReconCluster.hxx>
#include <TCaptLog.hxx>
#include <HEPUnits.hxx>
#include <TRuntimeParameters.hxx>

#include <memory>
#include <cmath>

namespace {
    struct HitProximity {
        double operator() (CP::THandle<CP::THit> lhs, 
                           CP::THandle<CP::THit> rhs) {
            return (lhs->GetPosition() - rhs->GetPosition()).Mag();
        }
    };

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
CP::TClusterSlice::BreakCluster(CP::THandle<CP::TReconBase> in) {
    CP::THandle<CP::TReconObjectContainer> result(
        new CP::TReconObjectContainer("result"));
    // This only looks at clusters.  If it's not a cluster, then just copy out
    // the input the input object.
    CP::THandle<CP::TReconCluster> cluster = in;
    if (!cluster) {
        result->push_back(in);
        return result;
    }

    /// Check the input hits.  If there are not enough, then just return the 
    CP::THandle<CP::THitSelection> inputHits = cluster->GetHits();
    if (!inputHits || inputHits->size() < fMinHits) {
        result->push_back(in);
        return result;
    }
    
    /// Copy the input hits and sort them by the Z position.
    std::auto_ptr<CP::THitSelection> hits(new THitSelection);
    hits->reserve(inputHits->size());
    std::copy(inputHits->begin(), inputHits->end(), std::back_inserter(*hits));
    std::sort(hits->begin(), hits->end(), HitZSort());
    
    // Distance between first and last point.
    double deltaZ = std::abs(hits->front()->GetPosition().Z() 
                             - hits->back()->GetPosition().Z());
    // The number of steps based the cluster step in Z.
    int steps = deltaZ/fClusterStep;
    // The adjusted step size so that each step is the same "thickness" in Z.
    double zStep = deltaZ / (steps+1);
    
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
        CP::THandle<CP::TReconCluster> cluster(new TReconCluster);
        // The hits between first and curr should be run through the density
        // cluster again since it's possible (not unlikely) for the track to
        // have curved back in Z.
        cluster->FillFromHits("zCluster", first,curr);
        result->push_back(cluster);
        // Reset first to start looking for a new set of hits.
        first = curr;
    }
    if (first != end) {
        // Build the final cluster.
        CP::THandle<CP::TReconCluster> cluster(new TReconCluster);
        // The hits between first and curr should be run through the density
        // cluster again since it's possible (not unlikely) for the track to
        // have curved back in Z.
        cluster->FillFromHits("zCluster", first,curr);
        result->push_back(cluster);
    }

    return result;
}

CP::THandle<CP::TReconObjectContainer>
CP::TClusterSlice::MakeTracks(CP::THandle<CP::TReconObjectContainer> input) {
    CP::THandle<CP::TReconObjectContainer> 
        output(new CP::TReconObjectContainer);
    CP::TReconObjectContainer remains;
    CP::TReconObjectContainer seed;
    std::copy(input->begin(), input->end(), std::back_inserter(remains));
    CP::TTrackFit trackFitter;
    for (int i=0; i<10; ++i) {
        CaptSevere("Look for seed");
        // Find a seed.
        seed.clear();
        std::auto_ptr<CP::TBestTube> bestTube(new CP::TBestTube);
        bestTube->Process(remains);
        bestTube->FillSeed(seed);
        bestTube->FillRemains(remains);
        if (seed.empty()) break;
        // Make the track.
        std::auto_ptr<CP::TLinearRoad> road(new CP::TLinearRoad);
        road->Process(seed,remains);
        CP::THandle<CP::TReconTrack> track = road->GetTrack();
        if (track) {
            // Fit the track and then save it to the output.
            track = trackFitter(track);
            output->push_back(track);
        }
        // Get the remaining hits.
        road->FillRemains(remains);
        // Check if we continue.
        if (remains.empty()) break;
    }
    return output;
}

CP::THandle<CP::TAlgorithmResult>
CP::TClusterSlice::Process(const CP::TAlgorithmResult& input,
                             const CP::TAlgorithmResult&,
                             const CP::TAlgorithmResult&) {

    CP::THandle<CP::TReconObjectContainer> inputObjects 
        = input.GetResultsContainer();
    if (!inputObjects) {
        CaptError("No input objects");
        return CP::THandle<CP::TAlgorithmResult>();
    }

    CaptLog("TClusterSlice Process " 
            << GetEvent().GetContext()
            << " w/ " <<  inputObjects->size() << " objects");


    CP::THandle<CP::TAlgorithmResult> result = CreateResult();
    std::auto_ptr<CP::TReconObjectContainer> 
        final(new CP::TReconObjectContainer("final"));
    std::auto_ptr<CP::THitSelection> used(new CP::THitSelection("used"));

    for (CP::TReconObjectContainer::iterator obj = inputObjects->begin();
         obj != inputObjects->end(); ++obj) {
        CP::THandle<CP::TReconObjectContainer> newObjs = BreakCluster(*obj);
        CP::THandle<CP::TReconObjectContainer> tracks = MakeTracks(newObjs);
        std::copy(tracks->begin(), tracks->end(), std::back_inserter(*final));
    }

    // Copy all of the hits that got added to a reconstruction object into the
    // used hit selection.
    CP::THandle<CP::THitSelection> hits 
        = CP::hits::ReconHits(final->begin(), final->end());
    if (hits) {
        used->reserve(hits->size());
        std::copy(hits->begin(), hits->end(), std::back_inserter(*used));
    }

    result->AddHitSelection(used.release());
    result->AddResultsContainer(final.release());

    return result;
}
