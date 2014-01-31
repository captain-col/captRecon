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

    /////////////////////////////////////////////////////////////////////
    /// Estimate the state at each node.
    /////////////////////////////////////////////////////////////////////
    // Hold the last cluster seen.
    CP::THandle<CP::TReconCluster> prevCluster;
    double energyDeposit = 0.0;
    for (TReconNodeContainer::iterator n = nodes.begin();
         n != nodes.end(); ++n) {
        // Find the next cluster in the track.
        CP::THandle<CP::TReconCluster> nextCluster;
        if ((n + 1) != nodes.end()) nextCluster = (*(n+1))->GetObject();

        // Get the state from the node.  It had better be a trackState!.
        CP::THandle<CP::TTrackState> trackState = (*n)->GetState();

        // Get the object from the node.  It had better be a cluster.
        CP::THandle<CP::TReconCluster> cluster = (*n)->GetObject();
        if (!cluster) {
            CaptError("Object is not a cluster");
            continue;
        }

        CP::THandle<CP::TReconCluster> first = prevCluster;
        if (!first) first = cluster;

        CP::THandle<CP::TReconCluster> last = nextCluster;
        if (!last) last = cluster;

        if (first == last) {
            CaptError("This 'cannot' happen. The track node has no neighbors.");
            CaptError("  Number of nodes " << nodes.size());
            CaptError("  At position " << n - nodes.begin());
            CaptError("  prevCluster " << prevCluster);
            CaptError("  first       " << first);
            CaptError("  cluster     " << cluster);
            CaptError("  last        " << last);
            CaptError("  nextCluster " << nextCluster);
            prevCluster = cluster;
            continue;
        }
        
        // Make a crude estimate of the direction and direction covariance.
        TVector3 dir = (last->GetPosition().Vect()-first->GetPosition().Vect());
        TMatrixD dirCov(3,3);
        double length2 = dir.Mag2();
        dir = dir.Unit();
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
        energyDeposit += cluster->GetEDeposit();
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
    trackFront->SetEDeposit(energyDeposit);

    // Set the back state of the track.
    CP::THandle<CP::TTrackState> trackBack = input->GetBack();
    CP::THandle<CP::TTrackState> lastNodeState = nodes.back()->GetState();
    *trackBack = *lastNodeState;
    trackBack->SetEDeposit(0.0);

    int trackDOF = 3*nodes.size() - 6;
    input->SetStatus(TReconBase::kSuccess);
    input->SetAlgorithmName("TSegmentTrackFit");
    input->SetQuality(0.0);
    input->SetNDOF(trackDOF);

    return input;
}
