#ifndef CreateTrack_hxx_seen
#define CreateTrack_hxx_seen

#include "TSegmentTrackFit.hxx"
#include "HitUtilities.hxx"
#include "ECaptRecon.hxx"

#include <TReconTrack.hxx>
#include <THandle.hxx>

namespace CP {

    /// A base exception for the create track template.
    EXCEPTION(ECreateTrack,ECaptRecon);

    /// An exception that there are repeated objects in the track.
    EXCEPTION(ETrackRepeatedObject, ECreateTrack);

    /// An exception that there are repeated objects in the track.
    EXCEPTION(ETrackNonCluster, ECreateTrack);

    /// Take iterators from a container holding TReconCluster objects in the
    /// right order for the track nodes, and construct a track.  The first
    /// argument becomes the name of the algorithm that created the track.
    /// The track will be "fitted" using the TSegmentTrackFit which actually
    /// just connects the clusters with straight line segments.  The resulting
    /// track can be refit by TTrackFit.
    template<typename clusterIterator>
    CP::THandle<CP::TReconTrack> 
    CreateTrack(const char* name, clusterIterator begin, clusterIterator end,
                bool verify=true);

    /// Take iterators from a container holding THandle<THit> objects (in no
    /// particular order), and create a track.  The first argument becomes the
    /// name of the algorithm that created the track.  The track will be
    /// "fitted" using the TSegmentTrackFit which actually just connects the
    /// clusters with straight line segments.  The resulting track can be
    /// refit by TTrackFit.
    template<typename hitIterator>
    CP::THandle<CP::TReconTrack> 
    CreateTrack(const char* name, hitIterator begin, hitIterator end,
                const TVector3& approxDir);

};

//////////////////////////////////////////////////////////////////
// IMPLEMENTATION
//////////////////////////////////////////////////////////////////

template<typename clusterIterator>
CP::THandle<CP::TReconTrack> 
CP::CreateTrack(const char* name, clusterIterator begin, clusterIterator end,
                bool verify=true) {

    if (verify) {
        for (clusterIterator i = begin; i!=end; ++i) {
            clusterIterator j = i;
            while ((++j) != end) {
                if (CP::GetPointer(*i) != CP::GetPointer(*j)) continue;
                CaptError("Invalid track: multiple copies of object");
                throw CP::ETrackRepeatedObject();
            }
        }
    }

    CP::THandle<CP::TReconTrack> track(new CP::TReconTrack);
    track->SetAlgorithmName(name);
    track->SetStatus(CP::TReconBase::kSuccess);
    track->AddDetector(CP::TReconBase::kTPC);
    track->SetName("track");

    std::set< CP::THandle<CP::THit> > hits;
    CP::hits::ReconHits(begin,end,hits);
    CP::THitSelection* trackHits = new CP::THitSelection("trackHits");
    std::copy(hits.begin(), hits.end(),std::back_inserter(*trackHits));
    track->AddHits(trackHits);

    TReconNodeContainer& nodes = track->GetNodes();
    while (begin != end) {
        CP::THandle<CP::TReconCluster> cluster = *begin;
        if (!cluster) {
            CaptError("Invalid track: object not a cluster");
            throw CP::ETrackNonCluster();
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

    // The segment fit was applied to make sure all of the fields are
    // initialized, but make the track as not fit anyway.
    track->ClearStatus(CP::TReconBase::kStatusMask);

    return track;
}

template<typename hitIterator>
CP::THandle<CP::TReconTrack> 
CP::CreateTrack(const char* name, hitIterator begin, hitIterator end,
                const TVector3& approxDir) {

    CP::THandle<CP::TReconObjectContainer> clusters
        = CreateClusters(name, begin, end, 
                         approxDir);

    return CreateTrack(name,clusters->begin(), clusters->end(), false);
}
#endif
