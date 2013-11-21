#include "TSegmentTrackFit.hxx"

#include <TReconNode.hxx>
#include <TReconCluster.hxx>
#include <HEPUnits.hxx>
#include <HEPConstants.hxx>
#include <TCaptLog.hxx>

#include <TMatrixD.h>
#include <TPrincipal.h>

#include <memory>

CP::TSegmentTrackFit::TSegmentTrackFit() {}
CP::TSegmentTrackFit::~TSegmentTrackFit() {}

CP::THandle<CP::TReconTrack>
CP::TSegmentTrackFit::Apply(CP::THandle<CP::TReconTrack>& input) {

    TReconNodeContainer& nodes = input->GetNodes();
    if (nodes.size() < 2) {
        CaptError("Not enough nodes to fit.");
        return CP::THandle<CP::TReconTrack>();
    }

    CaptLog("Start segment fit with " << nodes.size() << " nodes");

    /////////////////////////////////////////////////////////////////////
    /// Estimate the state at each node.
    /////////////////////////////////////////////////////////////////////
    // Hold the last cluster seen.
    CP::THandle<CP::TReconCluster> prevCluster;
    for (TReconNodeContainer::iterator n = nodes.begin();
         n != nodes.end(); ++n) {
        // Find the next cluster in the track.
        CP::THandle<CP::TReconCluster> nextCluster;
        if ((n + 1) != nodes.end()) nextCluster = (*(n+1))->GetState();

        // Get the state from the node.  It had better be a trackState!.
        CP::THandle<CP::TTrackState> trackState = (*n)->GetState();

        // Get the object from the node.  It had better be a cluster.
        CP::THandle<CP::TReconCluster> cluster = (*n)->GetObject();

        CP::THandle<CP::TReconCluster> first = prevCluster;
        if (!first) first = cluster;

        CP::THandle<CP::TReconCluster> last = nextCluster;
        if (!last) last = cluster;

        if (first == last) {
            CaptError("This 'cannot' happen. The track node has no neighbors.");
            continue;
        }
        
        // Make a crude estimate of the direction and direction covariance.
        TVector3 dir = (last->GetPosition().Vect()-first->GetPosition().Vect());
        TMatrixD dirCov(3,3);
        double length2 = dir.Mag2();
        CP::THandle<CP::TClusterState> firstState = first->GetState();
        CP::THandle<CP::TClusterState> lastState = last->GetState();
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                double v = firstState->GetPositionCovariance(i,j)
                    + lastState->GetPositionCovariance(i,j);
                v = v/length2;
                dirCov(i,j) = v;
            }
        }            

        // Set the track state.
        trackState->SetEDeposit(cluster->GetEDeposit());
        trackState->SetPosition(cluster->GetPosition().X(), 
                                cluster->GetPosition().Y(),
                                cluster->GetPosition().Z(),
                                cluster->GetPosition().T());
        trackState->SetDirection(dir.X(),dir.Y(),dir.Z());

        CP::THandle<CP::TClusterState> clusterState = cluster->GetState();
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                double v = clusterState->GetPositionCovariance(i,j);
                trackState->SetPositionCovariance(i,j,v);
                trackState->SetDirectionCovariance(i,j,dirCov(i,j));
            }
        }

        prevCluster = cluster;
    }

    // Set the front state of the track.
    CP::THandle<CP::TTrackState> trackFront = input->GetFront();
    CP::THandle<CP::TTrackState> firstNodeState = nodes.front()->GetState();
    *trackFront = *firstNodeState;

    // Now move the front state upstream to the position of the first hit.
    // Notice that the front state is not at the same location as the first
    // node.  
    if (nodes.front()->GetObject()) {
        // There is an object, and there had *better* be one or something is
        // going horribly wrong.
        TVector3 pos = trackFront->GetPosition().Vect();
        TVector3 dir = trackFront->GetDirection();
        CP::THandle<CP::THitSelection> hits 
            = nodes.front()->GetObject()->GetHits();
        if (!hits) {
            // There had also better be hits!
            CaptError("No hits in object!");
            abort();
        }
        double upstream = 0.0;
        for (CP::THitSelection::iterator h = hits->begin(); 
             h != hits->end(); ++h) {
            TVector3 diff = (*h)->GetPosition() - pos;
            upstream = std::min(upstream, diff*dir);
        }
        pos = pos + upstream*dir;
        trackFront->SetPosition(pos.X(), pos.Y(), pos.Z(),
                                trackFront->GetPosition().T());
    }

    // Set the back state of the track.
    CP::THandle<CP::TTrackState> trackBack = input->GetBack();
    CP::THandle<CP::TTrackState> lastNodeState = nodes.back()->GetState();
    *trackBack = *lastNodeState;

    // Now move the back state downstream to the position of the last hit.
    // See the comments for "trackFront".
    if (nodes.back()->GetObject()) {
        // There is an object, and there had *better* be one or something is
        // going horribly wrong.
        TVector3 pos = trackBack->GetPosition().Vect();
        TVector3 dir = trackBack->GetDirection();
        CP::THandle<CP::THitSelection> hits 
            = nodes.back()->GetObject()->GetHits();
        if (!hits) {
            // There had also better be hits!
            CaptError("No hits in object!");
            abort();
        }
        double downstream = 0.0;
        for (CP::THitSelection::iterator h = hits->begin(); 
             h != hits->end(); ++h) {
            TVector3 diff = (*h)->GetPosition() - pos;
            downstream = std::max(downstream, diff*dir);
        }
        pos = pos + downstream*dir;
        trackBack->SetPosition(pos.X(), pos.Y(), pos.Z(),
                               trackBack->GetPosition().T());
    }

    int trackDOF = 3*nodes.size() - 6;
    input->SetStatus(TReconBase::kSuccess);
    input->SetAlgorithmName("TSegmentTrackFit");
    input->SetQuality(0.0);
    input->SetNDOF(trackDOF);

    return input;
}
