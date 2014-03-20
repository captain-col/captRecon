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
    double energyVariance = 0.0;
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
#ifdef USE_POSITION_COVARIANCE_FOR_DIR_COVARIANCE
                double v = firstState->GetPositionCovariance(i,j)
                    + lastState->GetPositionCovariance(i,j);
#else
                double v = first->GetMoments()(i,j) + last->GetMoments()(i,j);
#endif
                v = v/length2;
                dirCov(i,j) = v;
            }
        }            

        // Set the track state.
        energyDeposit += cluster->GetEDeposit();
        energyVariance += cluster->GetEDepositVariance();
        trackState->SetEDeposit(cluster->GetEDeposit());
        trackState->SetEDepositVariance(cluster->GetEDepositVariance());
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
    CP::THandle<CP::TTrackState> frontState = input->GetFront();
    CP::THandle<CP::TTrackState> firstNodeState = nodes.front()->GetState();
    *frontState = *firstNodeState;
    frontState->SetEDeposit(energyDeposit);
    frontState->SetEDepositVariance(energyVariance);

    // Notice that the front state is not at the same location as the first
    // node.  Adjust the position of the front of the track to be at the
    // beginning of the first nodes object.
    TVector3 frontPos = frontState->GetPosition().Vect();
    TVector3 frontDir = frontState->GetDirection();
    CP::THandle<CP::THitSelection> frontHits 
        = nodes.front()->GetObject()->GetHits();
    if (!frontHits) {
        // There had also better be hits!
        CaptError("No hits in object!");
        abort();
    }
    double upstream = 0.0;
    for (CP::THitSelection::iterator h = frontHits->begin(); 
         h != frontHits->end(); ++h) {
        TVector3 diff = (*h)->GetPosition() - frontPos;
        upstream = std::min(upstream, diff*frontDir);
    }
    frontPos = frontPos + upstream*frontDir;
    frontState->SetPosition(frontPos.X(), frontPos.Y(), frontPos.Z(),
                            frontState->GetPosition().T());

    // Set the back state of the track.
    CP::THandle<CP::TTrackState> backState = input->GetBack();
    CP::THandle<CP::TTrackState> lastNodeState = nodes.back()->GetState();
    *backState = *lastNodeState;
    backState->SetEDeposit(0.0);
    backState->SetEDepositVariance(0.0);

    // Now move the back state downstream to the position of the last hit.
    // See the comments for "frontState".
    TVector3 backPos = backState->GetPosition().Vect();
    TVector3 backDir = backState->GetDirection();
    CP::THandle<CP::THitSelection> backHits 
        = nodes.back()->GetObject()->GetHits();
    if (!backHits) {
        // There had also better be hits!
        CaptError("No hits in object!");
        abort();
    }
    double downstream = 0.0;
    for (CP::THitSelection::iterator h = backHits->begin(); 
         h != backHits->end(); ++h) {
        TVector3 diff = (*h)->GetPosition() - backPos;
        downstream = std::max(downstream, diff*backDir);
    }
    backPos = backPos + downstream*backDir;
    backState->SetPosition(backPos.X(), backPos.Y(), backPos.Z(),
                           backState->GetPosition().T());

    // Finish up the fit information...
    int trackDOF = 3*nodes.size() - 6;
    input->SetStatus(TReconBase::kSuccess);
    input->SetStatus(TReconBase::kRan);
    input->SetAlgorithmName("TSegmentTrackFit");
    input->SetQuality(0.0);
    input->SetNDOF(trackDOF);

    return input;
}
