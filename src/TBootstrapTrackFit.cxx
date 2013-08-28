#include "TBootstrapTrackFit.hxx"

CP::TTrackFit::TTrackFit() {}
CP::TTrackFit::~TTrackFit() {}

CP::THandle<CP::TReconTrack>
CP::TTrackFit::Apply(CP::THandle<CP::TReconTrack>& input) {
    return input;
}
