#include "TBestTubeTrack.hxx"
#include "HitUtilities.hxx"
#include "TLinearRoad.hxx"
#include "TBestTube.hxx"
#include "TTrackFit.hxx"

#include <THandle.hxx>
#include <TReconCluster.hxx>
#include <TCaptLog.hxx>
#include <HEPUnits.hxx>
#include <TUnitsTable.hxx>
#include <TRuntimeParameters.hxx>

#include <memory>
#include <iostream>
#include <sstream>
#include <cmath>

CP::TBestTubeTrack::TBestTubeTrack()
    : TAlgorithm("TBestTubeTrack", "Find Tracks Based on Tubes") {
}

CP::TBestTubeTrack::~TBestTubeTrack() { }

CP::THandle<CP::TReconObjectContainer>
CP::TBestTubeTrack::MakeTracks(const CP::TReconObjectContainer& input,
                               CP::THandle<CP::TAlgorithmResult> thisResult) {
    CP::THandle<CP::TReconObjectContainer> 
        output(new CP::TReconObjectContainer);
    CP::TReconObjectContainer remains;
    CP::TReconObjectContainer seed;
    std::copy(input.begin(), input.end(), std::back_inserter(remains));
    CP::TTrackFit trackFitter;
    for (int i=0; i<50; ++i) {
        CaptNamedLog("TBestTubeTrack","Look for seed " << i
                     << " in " << remains.size() << " clusters");
        // Find a seed.
        seed.clear();
        std::auto_ptr<CP::TBestTube> bestTube(new CP::TBestTube);
        std::size_t lastRemains = remains.size();
        bestTube->Process(remains);
        bestTube->FillSeed(seed);
        bestTube->FillRemains(remains);
        CaptNamedLog("TBestTubeTrack","Found seed with " 
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
            CaptNamedLog("TBestTubeTrack","Track found with "
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
CP::TBestTubeTrack::Process(const CP::TAlgorithmResult& input,
                           const CP::TAlgorithmResult&,
                           const CP::TAlgorithmResult&) {

    CP::THandle<CP::TReconObjectContainer> inputObjects 
        = input.GetResultsContainer();

    if (!inputObjects) {
        CaptError("No input objects");
        return CP::THandle<CP::TAlgorithmResult>();
    }

    if (inputObjects->size() > 1000) {
        CaptError("Too many input objects");
        return CP::THandle<CP::TAlgorithmResult>();
    }

    CaptLog("TBestTubeTrack Process " << GetEvent().GetContext());
    
    // Create the output containers.
    CP::THandle<CP::TAlgorithmResult> result = CreateResult();
    std::auto_ptr<CP::TReconObjectContainer> 
        final(new CP::TReconObjectContainer("final"));
    std::auto_ptr<CP::THitSelection> used(new CP::THitSelection("used"));

    CP::THandle<CP::TReconObjectContainer> tracks 
        = MakeTracks(*inputObjects,result);
    CaptLog("   Objects found: " << tracks->size());
    std::copy(tracks->begin(),tracks->end(),std::back_inserter(*final));

    // Copy all of the hits that got added to a reconstruction object into the
    // used hit selection.
    CP::THandle<CP::THitSelection> hits 
        = CP::hits::ReconHits(final->begin(), final->end());
    if (hits) {
        used->reserve(hits->size());
        std::copy(hits->begin(), hits->end(), std::back_inserter(*used));
    }

    result->AddHits(used.release());
    result->AddResultsContainer(final.release());

    return result;
}
