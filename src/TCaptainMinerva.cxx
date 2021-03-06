#include "TCaptainMinerva.hxx"
#include "TCluster3D.hxx"
#include "TDensityCluster.hxx"
#include "TClusterSlice.hxx"
#include "TClusterMerge.hxx"
#include "TBestTubeTrack.hxx"
#include "TMinimalSpanningTrack.hxx"
#include "TSplitTracks.hxx"
#include "TMergeTracks.hxx"
#include "TDisassociateHits.hxx"
#include "TCombineOverlaps.hxx"
#include "TDestroyShortTracks.hxx"
#include "TClusterUnusedHits.hxx"

#include "HitUtilities.hxx"

#include <THandle.hxx>
#include <TReconHit.hxx>
#include <TReconCluster.hxx>
#include <TCaptLog.hxx>
#include <CaptGeomId.hxx>
#include <HEPUnits.hxx>

#include <algorithm>
#include <memory>
#include <cmath>

CP::TCaptainMinerva::TCaptainMinerva()
    : TAlgorithm("TCaptainMinerva", "Cluster Wire Hits") { }

CP::TCaptainMinerva::~TCaptainMinerva() { }

CP::THandle<CP::TAlgorithmResult>
CP::TCaptainMinerva::Process(const CP::TAlgorithmResult& driftInput,
                             const CP::TAlgorithmResult& pmtInput,
                             const CP::TAlgorithmResult&) {

    CaptLog("TCaptainMinerva Process " << GetEvent().GetContext());
    CP::THandle<CP::TAlgorithmResult> result(CP::TAlgorithm::CreateResult());

    ///////////////////////////////////////////////////////////////////
    // Get the  input hits.
    ///////////////////////////////////////////////////////////////////

    CP::THandle<CP::THitSelection> wires = driftInput.GetHits();
    if (!wires) {
        CaptError("No wire hits");
        return CP::THandle<CP::TAlgorithmResult>();
    }

    CP::THandle<CP::THitSelection> pmts = pmtInput.GetHits();

    // The final objects from this will be copied into the final recon
    // container.  It might be NULL if there is a problem.
    CP::THandle<CP::TAlgorithmResult> currentResult;

    // All the hits.  This is actually the result of the 3D clustering
    // algorithm.  This is set after the clustering, and might be NULL if
    // there is a problem.
    CP::THandle<CP::THitSelection> allHits;

    ///////////////////////////////////////////////////////////
    // Apply the sub-algorithms in sequence.
    ///////////////////////////////////////////////////////////
    do {
        // Find the time zero and the 3D hits.
        CP::THandle<CP::TAlgorithmResult> cluster3DResult;
        if (pmts) {
            cluster3DResult = Run<CP::TCluster3D>(*wires,*pmts);
        }
        else {
            cluster3DResult = Run<CP::TCluster3D>(*wires);
        }
        if (!cluster3DResult) break;
        currentResult = cluster3DResult;
        allHits = currentResult->GetHits();
        result->AddDatum(currentResult);

        // Not usually applied because it is replaced by TClusterSlice.

#ifdef Apply_TDensityCluster
        CP::THandle<CP::TAlgorithmResult> densityClusterResult
            = Run<CP::TDensityCluster>(*currentResult);
        if (!densityClusterResult) break;
        currentResult = densityClusterResult;
        result->AddDatum(currentResult);
#endif
	
#define Apply_TClusterSlice
#ifdef Apply_TClusterSlice
        CP::THandle<CP::TAlgorithmResult> clusterSliceResult
            = Run<CP::TClusterSlice>(*currentResult);
        if (!clusterSliceResult) break;
        currentResult = clusterSliceResult;
        result->AddDatum(currentResult);
#endif

#define Apply_TClusterMerge
#ifdef Apply_TClusterMerge
        CP::THandle<CP::TAlgorithmResult> clusterMergeResult
            = Run<CP::TClusterMerge>(*currentResult);
        if (!clusterMergeResult) break;
        currentResult = clusterMergeResult;
        result->AddDatum(currentResult);
#endif


#define Apply_TMinimalSpanningTrack
#ifdef Apply_TMinimalSpanningTrack
        CP::THandle<CP::TAlgorithmResult> minimalSpanningResult
            = Run<CP::TMinimalSpanningTrack>(*currentResult);
        if (!minimalSpanningResult) break;
        currentResult = minimalSpanningResult;
        result->AddDatum(currentResult);
#endif

#define Apply_TSplitTracks
#ifdef Apply_TSplitTracks
        CP::THandle<CP::TAlgorithmResult> splitTracksResult
            = Run<CP::TSplitTracks>(*currentResult);
        if (!splitTracksResult) break;
        currentResult = splitTracksResult;
        result->AddDatum(currentResult);
#endif
	
#define Apply_TMergeTracks
#ifdef Apply_TMergeTracks
        CP::THandle<CP::TAlgorithmResult> mergeTracksResult
            = Run<CP::TMergeTracks>(*currentResult);
        if (!mergeTracksResult) break;
        currentResult = mergeTracksResult;
        result->AddDatum(currentResult);
#endif

#define Apply_TDisassociateHits
#ifdef Apply_TDisassociateHits
        CP::THandle<CP::TAlgorithmResult> disassociateHitsResult
            = Run<CP::TDisassociateHits>(*currentResult);
        if (!disassociateHitsResult) break;
        currentResult = disassociateHitsResult;
        result->AddDatum(currentResult);
#endif

#define Apply_TCombineOverlaps
#ifdef Apply_TCombineOverlaps
        CP::THandle<CP::TAlgorithmResult> combineOverlapsResult
            = Run<CP::TCombineOverlaps>(*currentResult);
        if (!combineOverlapsResult) break;
        currentResult = combineOverlapsResult;
        result->AddDatum(currentResult);
#endif

#define Destroy_Short_Tracks
#ifdef Destroy_Short_Tracks
        //Destroy tracks with reconstructed length l < 15mm
        CP::THandle<CP::TAlgorithmResult> destroyShortTracks
            = Run<CP::TDestroyShortTracks>(*currentResult);
        if (!destroyShortTracks) break;
        currentResult = destroyShortTracks;
        result->AddDatum(currentResult);	
#endif

#define Cluster_Unused_Hits
#ifdef Cluster_Unused_Hits
        // Cluster unused by this time 3D hits by position. 
        CP::THandle<CP::TAlgorithmResult> clusterUnusedHits
            = Run<CP::TClusterUnusedHits>(*currentResult,*allHits);
        if (!clusterUnusedHits) break;
        currentResult = clusterUnusedHits;
        result->AddDatum(currentResult);	
#endif
	
    } while (false);
    
    if (!currentResult) {
        CaptError("Reconstruct failed");
        return result;
    }

    std::unique_ptr<CP::TReconObjectContainer> 
        finalObjects(new CP::TReconObjectContainer("final"));

    // Copy the reconstruction objects from the last algorithm result into the
    // final reconstruction objects for this result.
    CP::THandle<CP::TReconObjectContainer> currentObjects
        = currentResult->GetResultsContainer("final");
    if (currentObjects) {
        std::copy(currentObjects->begin(), currentObjects->end(),
                  std::back_inserter(*finalObjects));
    }
    // Save the hits, but only if there aren't to many.
    if (allHits && allHits->size() < 10000) {

     

        std::unique_ptr<CP::THitSelection> 
            used(new CP::THitSelection("used"));
        std::unique_ptr<CP::THitSelection> 
            unused(new CP::THitSelection("unused"));
        // Get all of the hits in the final object and add them to used.
        CaptLog("Fill the used hits");
        CP::THandle<CP::THitSelection> hits 
            = CP::hits::ReconHits(finalObjects->begin(), finalObjects->end());
        if (hits) {
            std::copy(hits->begin(), hits->end(), std::back_inserter(*used));
        }
        
        CaptLog("Fill the unused hits");
        unused->reserve(allHits->size());
        std::copy(allHits->begin(), allHits->end(), 
                  std::back_inserter(*unused));
        CP::hits::Subtract(*unused,*used);
        result->AddHits(unused.release());
        result->AddHits(used.release());
 
    }

    CaptLog("Save the results");
    result->AddResultsContainer(finalObjects.release());
    //CP::THandle<CP::TReconObjectContainer> m3Objects
    //  = currentResult->GetResultsContainer("match3Tr");
    // result->AddResultsContainer(&(*m3Objects));

    return result;
}
