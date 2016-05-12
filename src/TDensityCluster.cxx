#include "TDensityCluster.hxx"
#include "HitUtilities.hxx"
#include "TPositionDensityCluster.hxx"
#include "CreateCluster.hxx"

#include <THandle.hxx>
#include <TReconCluster.hxx>
#include <TCaptLog.hxx>
#include <HEPUnits.hxx>
#include <TRuntimeParameters.hxx>

#include <memory>
#include <cmath>

CP::TDensityCluster::TDensityCluster()
    : TAlgorithm("TDensityCluster", "Find Simply Connected Hits") {

    fMinPoints = CP::TRuntimeParameters::Get().GetParameterI(
        "captRecon.densityCluster.minPoints");

    fMaxDist = CP::TRuntimeParameters::Get().GetParameterD(
        "captRecon.densityCluster.maxDistance");


}

CP::TDensityCluster::~TDensityCluster() { }

CP::THandle<CP::TAlgorithmResult>
CP::TDensityCluster::Process(const CP::TAlgorithmResult& input,
                             const CP::TAlgorithmResult&,
                             const CP::TAlgorithmResult&) {

    CP::THandle<CP::THitSelection> inputHits = input.GetHits();
    if (!inputHits) {
        CaptError("No input hits");
        return CP::THandle<CP::TAlgorithmResult>();
    }

    CaptLog("TDensityCluster Process " 
            << GetEvent().GetContext()
            << " w/ " <<  inputHits->size() << " hits");


    CP::THandle<CP::TAlgorithmResult> result = CreateResult();
    std::unique_ptr<CP::TReconObjectContainer> 
        final(new CP::TReconObjectContainer("final"));
    std::unique_ptr<CP::THitSelection> used(new CP::THitSelection("used"));

    // This is the (dramatically) faster clustering algorithm for hits.
    typedef CP::TPositionDensityCluster< CP::THandle<CP::THit> > 
        ClusterAlgorithm;

    std::unique_ptr<ClusterAlgorithm> 
        clusterAlgorithm(new ClusterAlgorithm(fMinPoints,fMaxDist));

    clusterAlgorithm->Cluster(inputHits->begin(), inputHits->end());

    int nClusters = clusterAlgorithm->GetClusterCount();
    for (int i=0; i<nClusters; ++i) {
        const ClusterAlgorithm::Points& points 
            = clusterAlgorithm->GetCluster(i);
        CP::THandle<CP::TReconCluster> cluster
            = CreateCluster("TDensityCluster",points.begin(),points.end());
        CaptLog("   Cluster with " << cluster->GetHits()->size() << " hits");
        final->push_back(cluster);
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
    result->AddResultsContainer(final.release());

    return result;
}
