#include "TClusterSlice.hxx"
#include "HitUtilities.hxx"
#include "THitProximityCluster.hxx"

#include <THandle.hxx>
#include <TReconCluster.hxx>
#include <TCaptLog.hxx>
#include <HEPUnits.hxx>
#include <TUnitsTable.hxx>
#include <TRuntimeParameters.hxx>

#include <memory>
#include <iostream>
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

    CaptNamedLog("TClusterSlice","Make a maximum of " << steps << " slices"
                 << " in " 
                 << unit::AsString(deltaZ,"length") << " mm"
                 << " (" << unit::AsString(zStep,"length") << " mm each)");

    std::auto_ptr<CP::HitProximity::Cluster> 
        clusterAlgorithm(new CP::HitProximity::Cluster(1,6*unit::mm));

    int trials = 0;
    CP::THitSelection::iterator curr = hits->begin();
    CP::THitSelection::iterator end = hits->end();
    CP::THitSelection::iterator first = curr;
    while (curr != end) {
        // Make sure each cluster has at least 2*fMinPoints)
        if (curr-first < 2*fMinPoints) {
            ++curr;
            continue;
        }

        deltaZ = std::abs((*curr)->GetPosition().Z()
                          -(*first)->GetPosition().Z());
        // Make sure the cluster covers a minimum range in Z (determined by
        // hit resolution).
        double minZStep = 0.75*unit::mm;
        // Make sure that the cluster covers a range in Z.
        if (deltaZ < minZStep) {
            ++curr;
            continue;
        }

        // Make sure that the cluster at most the expected cluster step, but
        // limit the number of hits.
        if (deltaZ < zStep && curr-first < 200) {
            ++curr;
            continue;
        }

        ++trials;
        if (!(trials % 20)) {
            CaptError("Working on slice " << trials);
        }

        // Time for a new slice of clusters.
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
        CaptNamedInfo("TClusterSlice",
                     trials
                     << " -- Slice with " << nClusters
                     << " clusters from " << curr-first << " in slice"
                     << "   dZ: " << deltaZ);
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
            CaptLog("   Clusters made: " << cuts->size());
            std::copy(cuts->begin(), cuts->end(), std::back_inserter(*final));
        }
    }
    else if (inputHits) {
        CaptLog("   Using " <<  inputHits->size() << " hits");
        CP::THandle<CP::TReconObjectContainer> cuts = MakeSlices(inputHits);
        CaptLog("   Clusters made: " << cuts->size());
        std::copy(cuts->begin(), cuts->end(), std::back_inserter(*final));
    }


    // Copy all of the hits that got added to a reconstruction object into the
    // used hit selection.  But only if there aren't to many input hits.
    if (inputHits->size() < 5000) {
        CP::THandle<CP::THitSelection> hits 
            = CP::hits::ReconHits(final->begin(), final->end());
        std::auto_ptr<CP::THitSelection> used(new CP::THitSelection("used"));
        if (hits) {
            used->reserve(hits->size());
            std::copy(hits->begin(), hits->end(), std::back_inserter(*used));
        }
        result->AddHits(used.release());
    }

    result->AddResultsContainer(final.release());

    return result;
}
