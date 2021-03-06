#include "TClusterSlice.hxx"
#include "HitUtilities.hxx"
#include "TPositionDensityCluster.hxx"
#include "CreateCluster.hxx"
#include "ApproxArgonProperties.hxx"

#include <THandle.hxx>
#include <TReconCluster.hxx>
#include <TReconHit.hxx>
#include <TCaptLog.hxx>
#include <HEPUnits.hxx>
#include <TUnitsTable.hxx>
#include <TRuntimeParameters.hxx>

#include <memory>
#include <set>
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
    fMinHits = CP::TRuntimeParameters::Get().GetParameterD(
        "captRecon.clusterSlice.minHits");
    fMinStep = CP::TRuntimeParameters::Get().GetParameterD(
        "captRecon.clusterSlice.minStep");
    fClusterStep = CP::TRuntimeParameters::Get().GetParameterD(
        "captRecon.clusterSlice.clusterStep");

    fClusterExtent = CP::TRuntimeParameters::Get().GetParameterD(
        "captRecon.clusterSlice.clusterExtent");
    fClusterGrowth = CP::TRuntimeParameters::Get().GetParameterD(
        "captRecon.clusterSlice.clusterGrowth");
    fClusterCharge = 0.5*unit::mm*approxArgon::dEdX*approxArgon::Electrons;
}

CP::TClusterSlice::~TClusterSlice() { }

CP::THandle<CP::TReconObjectContainer> 
CP::TClusterSlice::MakeSlices(CP::THandle<CP::THitSelection> inputHits) {
    CP::THandle<CP::TReconObjectContainer> result(
        new CP::TReconObjectContainer("result"));

    if (!inputHits) return result;
    if (inputHits->size() < fMinHits) return result;
    
    /// Copy the input hits and sort them by the Z position.
    std::unique_ptr<CP::THitSelection> hits(new THitSelection);
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

    CaptNamedInfo("TClusterSlice","Make a maximum of " << steps << " slices"
                 << " in " 
                 << unit::AsString(deltaZ,"length") 
                 << " (" << unit::AsString(zStep,"length") << " each)");

    // This is the (dramatically) faster clustering algorithm for hits.
    typedef CP::TPositionDensityCluster<CP::THandle<CP::THit> > 
        ClusterAlgorithm;

// #define USE_SLICE_MINSIZE
#ifdef USE_SLICE_MINSIZE
    double eventScale = std::log(1.0*hits->size()+1.0)/std::log(fClusterGrowth);
    int minSize = eventScale;
    minSize = std::max(1,minSize);
#else
    int minSize = 1;
#endif
    
    int trials = 0;
    CP::THitSelection::iterator curr = hits->begin();
    CP::THitSelection::iterator end = hits->end();
    CP::THitSelection::iterator first = curr;
    while (curr != end) {
        // Make sure each slice has at least fMinPoints)
        if (curr-first < fMinPoints) {
            ++curr;
            continue;
        }

        // Check the case where the minimum cluster size (minSize) is bigger
        // than the minimum required points in a layer.
        if (curr-first <= minSize) {
            ++curr;
            continue;
        }

        deltaZ = std::abs((*curr)->GetPosition().Z()
                          -(*first)->GetPosition().Z());

        // Make sure that the cluster covers a range in Z.
        if (deltaZ < fMinStep) {
            ++curr;
            continue;
        }

        // Make sure that the cluster at most the expected cluster step, but
        // limit the number of hits.
        if (deltaZ < zStep && curr-first < 5000) {
            ++curr;
            continue;
        }

        ++trials;
        if (!(trials % 20)) {
            CaptNamedInfo("TClusterSlice", "Working on slice " << trials);
        }

        // Time for a new slice of clusters.
        ++curr;

        // Check that there are enough hits left for a new cluster.  If not,
        // then add all the remaining points to this cluster. (not really a
        // good plan, but it's got to suffice for now...}
        if (end-curr < fMinPoints) break;
        if (end-curr <= minSize) break;

        // The hits between first and curr should be run through the density
        // cluster again since it's very likely that the hits are disjoint in
        // the Z slice.
        std::unique_ptr<ClusterAlgorithm> 
            clusterAlgorithm(new ClusterAlgorithm(minSize,
                                                  fClusterExtent));
        clusterAlgorithm->Cluster(first,curr);
        int nClusters = clusterAlgorithm->GetClusterCount();
        CaptNamedInfo("TClusterSlice",
                     trials
                     << " -- Slice with " << nClusters
                     << " clusters from " << curr-first << " hits in slice"
                     << "   dZ: " << deltaZ);
        for (int i=0; i<nClusters; ++i) {
            const ClusterAlgorithm::Points& points 
                = clusterAlgorithm->GetCluster(i);
            CaptNamedVerbose("TClusterSlice","       Cluster " << i
                         << " with " << points.size() << " hits");
            CP::THandle<CP::TReconCluster> cluster
                = CreateCluster("zCluster",points.begin(),points.end());
            if (!cluster) continue;
            if (cluster->GetEDeposit() < fClusterCharge) {
                continue;
            }
            result->push_back(cluster);
        }
        // Reset first to start looking for a new set of hits.
        first = curr;
    }
    if (first != end) {
        // Build the final clusters.
        std::unique_ptr<ClusterAlgorithm> 
            clusterAlgorithm(new ClusterAlgorithm((int) minSize,
                                                  fClusterExtent));
        clusterAlgorithm->Cluster(first,end);
        int nClusters = clusterAlgorithm->GetClusterCount();
        CaptNamedVerbose("TClusterSlice","    Final slice with " << nClusters
                     << " clusters from " << end-first << " hits");
        for (int i=0; i<nClusters; ++i) {
            const ClusterAlgorithm::Points& points 
                = clusterAlgorithm->GetCluster(i);
            CaptNamedVerbose("TClusterSlice","       Cluster " << i
                         << " with " << points.size() << " hits");
            CP::THandle<CP::TReconCluster> cluster
                = CreateCluster("zCluster",points.begin(),points.end());
            if (!cluster) continue;
            if (cluster->GetEDeposit() < fClusterCharge) {
                continue;
            }
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

#ifdef CHECK_FOR_DUPLICATE_HITS
    CP::THandle<CP::THit> last;
    for (CP::THitSelection::iterator h = inputHits->begin();
         h != inputHits->end(); ++h) {
        do {
            if (!last) break;
            if (last->GetConstituentCount()!=(*h)->GetConstituentCount()) break;
            int match = 0;
            for (int i=0; i<last->GetConstituentCount(); ++i) {
                if (last->GetConstituent(i)->GetChannelId()
                    == (*h)->GetConstituent(i)->GetChannelId()) ++match;
            }
            if (match != last->GetConstituentCount()) break;
            if (std::abs(last->GetPosition().Z()-(*h)->GetPosition().Z())>1)
                break;
            CaptError("Duplicate hit " << match);
            for (int i=0; i<last->GetConstituentCount(); ++i) {
                CaptError("  Channels: "
                          << last->GetConstituent(i)->GetChannelId());
            }
            CaptError("  hit Z " << (*h)->GetPosition().Z());
            CaptError("  last Z " << last->GetPosition().Z());
        }  while (false);
        last = *h;
    }
#endif
    
    CaptLog("TClusterSlice Process " << GetEvent().GetContext());
    
    // Create the output containers.
    CP::THandle<CP::TAlgorithmResult> result = CreateResult();
    std::unique_ptr<CP::TReconObjectContainer> 
        final(new CP::TReconObjectContainer("final"));

    if (inputObjects) {
        for (CP::TReconObjectContainer::iterator obj = inputObjects->begin();
             obj != inputObjects->end(); ++obj) {
            CP::THandle<CP::THitSelection> hits = (*obj)->GetHits(); 
            if (!hits) {
                CaptError("Input object without hits");
                continue;
            }
            CP::THandle<CP::TReconObjectContainer> cuts = MakeSlices(hits);
            CaptLog("   Clusters made: " << cuts->size());
            std::copy(cuts->begin(), cuts->end(), std::back_inserter(*final));
        }
    }
    else if (inputHits) {
        CP::THandle<CP::TReconObjectContainer> cuts = MakeSlices(inputHits);
        CaptLog("   Clusters made: " << cuts->size());
        std::copy(cuts->begin(), cuts->end(), std::back_inserter(*final));
    }


    // Copy all of the hits that got added to a reconstruction object into the
    // used hit selection.  But only if there aren't to many input hits.
    if (inputHits->size() < 5000) {
        CP::THandle<CP::THitSelection> hits 
            = CP::hits::ReconHits(final->begin(), final->end());
        std::unique_ptr<CP::THitSelection> used(new CP::THitSelection("used"));
        if (hits) {
            used->reserve(hits->size());
            std::copy(hits->begin(), hits->end(), std::back_inserter(*used));
        }
        result->AddHits(used.release());
    }

    result->AddResultsContainer(final.release());

    return result;
}
