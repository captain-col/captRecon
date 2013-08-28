#include "TTrackFit.hxx"
#include "TBootstrapTrackFit.hxx"


CP::TTrackFit::TTrackFit()
    : fBootstrap(NULL) {}
CP::TTrackFit::~TTrackFit() {}

CP::THandle<CP::TReconTrack>
CP::TTrackFit::Apply(CP::THandle<CP::TReconTrack>& input) {
    // For now, the bootstrap fitter is the only fitter being used.
    if (!fBootstrap) {
        fBootstrap = new TBootstrapTrackFit;
    }
    return fBootstrap->Apply(input);
}
