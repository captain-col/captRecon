#include "TTrackFit.hxx"
#include "TBootstrapTrackFit.hxx"
#include "TClusterTrackFit.hxx"

CP::TTrackFit::TTrackFit()
    : fBootstrap(NULL), fCluster(NULL) {}
CP::TTrackFit::~TTrackFit() {}

CP::THandle<CP::TReconTrack>
CP::TTrackFit::Apply(CP::THandle<CP::TReconTrack>& input) {
    CP::THandle<CP::TReconTrack> result;
    if (input->GetNodes().size() > 10) {
        // For now, the bootstrap fitter is the only fitter being used.
        if (!fBootstrap) {
            fBootstrap = new TBootstrapTrackFit;
        }
        result = fBootstrap->Apply(input);
    }

    // If the bootstrap was successful, return it.
    if (result) return result;

    // The backstop fitter (the one for anything that can't be fit otherwise)
    // is the TClusterTrackFit.  That fit is super robust, but very
    // inaccurate.
    if (!fCluster) {
        fCluster = new TClusterTrackFit;
    }
    return fCluster->Apply(input);
}
