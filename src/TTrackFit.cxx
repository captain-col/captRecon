#include "TTrackFit.hxx"
#include "TBootstrapTrackFit.hxx"
#include "TClusterTrackFit.hxx"
#include "TSegmentTrackFit.hxx"

CP::TTrackFit::TTrackFit(int bootstrapIterations)
    : fBootstrap(NULL), fCluster(NULL), fSegment(NULL),
      fBootstrapIterations(bootstrapIterations) {}
CP::TTrackFit::~TTrackFit() {}

CP::THandle<CP::TReconTrack>
CP::TTrackFit::Apply(CP::THandle<CP::TReconTrack>& input) {
    CP::THandle<CP::TReconTrack> result;
    if (input->GetNodes().size() > 5) {
        // For longer tracks, the bootstrap fitter makes the best use of the
        // track information and produces a good result that follows the
        // multiple scattering of the track.
        if (!fBootstrap) {
            if (fBootstrapIterations>0) {
                fBootstrap = new TBootstrapTrackFit(fBootstrapIterations);
            }
            else {
                fBootstrap = new TBootstrapTrackFit;
            }
        }
        result = fBootstrap->Apply(input);
    }

    // If the bootstrap was successful, return it.
    if (result) return result;

    if (input->GetNodes().size() > 2) {
        // For shorter tracks, or if the boot strap fails, use the cluster
        // fitter.  This fits a straight line to the track and is super
        // robust, but not very accurate.
        if (!fCluster) {
            fCluster = new TClusterTrackFit;
        }
        result = fCluster->Apply(input);
    }

    // If the cluster fit was successful, return it.
    if (result) return result;

    if (!fSegment) {
        fSegment = new TSegmentTrackFit;
    }

    result = fSegment->Apply(input);

    return result;
}
