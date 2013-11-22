#ifndef CreateTrack_hxx_seen
#define CreateTrack_hxx_seen

#include "TSegmentTrackFit.hxx"
#include "HitUtilities.hxx"

#include <TReconTrack.hxx>
#include <THandle.hxx>

namespace CP {

    /// Take iterators from a container holding TReconCluster objects in the
    /// right order for the track nodes, and construct a track.  The first
    /// argument becomes the name of the algorithm that created the track.
    /// The track will be "fitted" using the TSegmentTrackFit which actually
    /// just connects the clusters with straight line segments.  The resulting
    /// track can be refit by TTrackFit.
    template<typename iterator>
    CP::THandle<CP::TReconTrack> 
    CreateTrack(const char* name, iterator begin, iterator end) {
        CP::THandle<CP::TReconTrack> track(new CP::TReconTrack);
        track->SetAlgorithmName(name);
        track->SetStatus(CP::TReconBase::kSuccess);
        track->AddDetector(CP::TReconBase::kTPC);
        track->SetName("track");

        CP::THandle<CP::THitSelection> hits 
            = CP::hits::ReconHits(begin,end);
        CP::THitSelection* trackHits = new CP::THitSelection("trackHits");
        std::copy(hits->begin(), hits->end(),std::back_inserter(*trackHits));
        track->AddHits(trackHits);

        TReconNodeContainer& nodes = track->GetNodes();
        while (begin != end) {
            CP::THandle<CP::TReconCluster> cluster = *begin;
            if (!cluster) {
                CaptError("Create track must have iterators to clusters");
                return CP::THandle<CP::TReconTrack>();
            }
            CP::THandle<CP::TReconNode> node(new CP::TReconNode);
            CP::THandle<CP::TReconState> state(new CP::TTrackState);
            CP::THandle<CP::TReconBase> object = cluster;
            node->SetState(state);
            node->SetObject(object);
            nodes.push_back(node);
            ++begin;
        }

        if (nodes.size() < 2) {
            CaptError("Not enough clusters for a track");
            return CP::THandle<CP::TReconTrack>();
        }

        CP::TSegmentTrackFit fitter;
        track = fitter.Apply(track);

        return track;
    }
};

#endif
