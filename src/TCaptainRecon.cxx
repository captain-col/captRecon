#include "TCaptainRecon.hxx"
#include "TCluster3D.hxx"
#include "TDensityCluster.hxx"
#include "TClusterSlice.hxx"
#include "TBestTubeTrack.hxx"

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

CP::TCaptainRecon::TCaptainRecon()
    : TAlgorithm("TCaptainRecon", "Cluster Wire Hits") { }

CP::TCaptainRecon::~TCaptainRecon() { }

CP::THandle<CP::TAlgorithmResult>
CP::TCaptainRecon::Process(const CP::TAlgorithmResult& driftInput,
                           const CP::TAlgorithmResult& pmtInput,
                           const CP::TAlgorithmResult&) {

    CaptLog("TCaptainRecon Process " << GetEvent().GetContext());
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
    if (!pmts) {
        CaptError("No PMT hits");
        return CP::THandle<CP::TAlgorithmResult>();
    }

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
        CP::THandle<CP::TAlgorithmResult> cluster3DResult
            = Run<CP::TCluster3D>(*wires,*pmts);
        if (!cluster3DResult) break;
        currentResult = cluster3DResult;
        allHits = currentResult->GetHits();
        result->AddDatum(currentResult);

#ifdef Apply_TDensityCluster
        // Cluster the 3D hits by position to find object candidates. 
        CP::THandle<CP::TAlgorithmResult> densityClusterResult
            = Run<CP::TDensityCluster>(*currentResult);
        if (!densityClusterResult) break;
        currentResult = densityClusterResult;
        result->AddDatum(currentResult);
#endif

#define Apply_TClusterSlice
#ifdef Apply_TClusterSlice
        // Cluster the 3D hits by position to find object candidates. 
        CP::THandle<CP::TAlgorithmResult> clusterSliceResult
            = Run<CP::TClusterSlice>(*currentResult);
        if (!clusterSliceResult) break;
        currentResult = clusterSliceResult;
        result->AddDatum(currentResult);
#endif

#define Apply_TBestTubeTrack
#ifdef Apply_TBestTubeTrack
        // Cluster the 3D hits by position to find object candidates. 
        CP::THandle<CP::TAlgorithmResult> bestTubeResult
            = Run<CP::TBestTubeTrack>(*currentResult);
        if (!bestTubeResult) break;
        currentResult = bestTubeResult;
        result->AddDatum(currentResult);
#endif

    } while (false);
    
    if (!currentResult) {
        CaptError("Reconstruct failed");
        return result;
    }
    
    std::auto_ptr<CP::TReconObjectContainer> 
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
        std::auto_ptr<CP::THitSelection> 
            used(new CP::THitSelection("used"));
        std::auto_ptr<CP::THitSelection> 
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

    return result;
}
