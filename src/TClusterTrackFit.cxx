#include "TClusterTrackFit.hxx"

#include <TReconNode.hxx>
#include <TReconCluster.hxx>
#include <HEPUnits.hxx>
#include <HEPConstants.hxx>
#include <TCaptLog.hxx>

// #include <ostreamTVector3.hxx>
// #include <ostreamTLorentzVector.hxx>

#include <TMatrixD.h>
#include <TPrincipal.h>

#include <memory>

namespace {

    // Take a position and turn it into a principal component value.
    double FindPrincipal(TPrincipal& pca, const TVector3& position) {
        double X[3] = {position.X(), position.Y(), position.Z()};
        double P[3] = {0,0,0};
        pca.X2P(X,P);
        return P[0];
    }

    // Take a principal component value and turn it into a position.
    TVector3 
    FindPosition(TPrincipal& pca, double principal) {
        double X[3] = {0,0,0};
        double P[3] = {principal,0,0};
        pca.P2X(P,X,3);
        return TVector3(X[0],X[1],X[2]);
    }
}

CP::TClusterTrackFit::TClusterTrackFit() {}
CP::TClusterTrackFit::~TClusterTrackFit() {}

CP::THandle<CP::TReconTrack>
CP::TClusterTrackFit::Apply(CP::THandle<CP::TReconTrack>& input) {

    TReconNodeContainer& nodes = input->GetNodes();
    if (nodes.size() < 2) {
        CaptError("Not enough nodes to fit.");
        return CP::THandle<CP::TReconTrack>();
    }

    CaptLog("Start cluster fit with " << nodes.size() << " nodes");

    /////////////////////////////////////////////////////////////////////
    /// Find the center position and direction of the track from the nodes.
    /////////////////////////////////////////////////////////////////////
    std::auto_ptr<TPrincipal> pca(new TPrincipal(3,""));
    for (TReconNodeContainer::iterator n = nodes.begin();
         n != nodes.end(); ++n) {
        // Get the object from the node.  It had better be a cluster.
        CP::THandle<CP::TReconCluster> cluster = (*n)->GetObject();
        if (!cluster) {
            CaptError("Missing cluster in track");
            return CP::THandle<CP::TReconTrack>();
        }
        double row1[3] 
            = {cluster->GetPosition().X()-cluster->GetLongAxis().X(),
               cluster->GetPosition().Y()-cluster->GetLongAxis().Y(),
               cluster->GetPosition().Z()-cluster->GetLongAxis().Z()};
        double row2[3] 
            = {cluster->GetPosition().X(),
               cluster->GetPosition().Y(),
               cluster->GetPosition().Z()};
        double row3[3] 
            = {cluster->GetPosition().X()+cluster->GetLongAxis().X(),
               cluster->GetPosition().Y()+cluster->GetLongAxis().Y(),
               cluster->GetPosition().Z()+cluster->GetLongAxis().Z()};
        for (double q = cluster->GetEDeposit(); q > 0; q -= 1000) {
            pca->AddRow(row1);
            pca->AddRow(row2);
            pca->AddRow(row2);
            pca->AddRow(row2);
            pca->AddRow(row3);
        }
    }
    pca->MakePrincipals();

    /// Find the front and last positions of the track and then find the track
    /// direction.
    CP::THandle<CP::TReconCluster> frontCluster = nodes.front()->GetObject();
    double frontPrincipal 
        = FindPrincipal(*pca, frontCluster->GetPosition().Vect());
    TVector3 frontPosition = FindPosition(*pca,frontPrincipal);
    
    CP::THandle<CP::TReconCluster> backCluster = nodes.back()->GetObject();
    double backPrincipal 
        = FindPrincipal(*pca, backCluster->GetPosition().Vect());
    TVector3 backPosition = FindPosition(*pca,backPrincipal);

    // The track direction is just the direction between the front and back
    // ends of the track.
    TVector3 dir = (backPosition-frontPosition).Unit();
    TMatrixD dirCov(3,3);
    double length2 = (backPosition-frontPosition).Mag2();
    CP::THandle<CP::TClusterState> frontState = frontCluster->GetState();
    CP::THandle<CP::TClusterState> backState = backCluster->GetState();
    for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
            double v = frontState->GetPositionCovariance(i,j)
                + backState->GetPositionCovariance(i,j);
            v = v/length2;
            dirCov(i,j) = v;
        }
    }            

    /////////////////////////////////////////////////////////////////////
    /// Estimate the state at each node.
    /////////////////////////////////////////////////////////////////////
    double logLikelihood = 0.0;
    for (TReconNodeContainer::iterator n = nodes.begin();
         n != nodes.end(); ++n) {
        // Get the state from the node.  It had better be a trackState!.
        CP::THandle<CP::TTrackState> trackState = (*n)->GetState();

        // Get the object from the node.  It had better be a cluster.
        CP::THandle<CP::TReconCluster> cluster = (*n)->GetObject();

        double nodePrincipal 
            = FindPrincipal(*pca, cluster->GetPosition().Vect());
        TVector3 nodePosition = FindPosition(*pca,nodePrincipal);

        // Set the track state.
        trackState->SetEDeposit(cluster->GetEDeposit());
        trackState->SetPosition(nodePosition.X(), 
                                nodePosition.Y(),
                                nodePosition.Z(),
                                cluster->GetPosition().T());
        trackState->SetDirection(dir.X(),dir.Y(),dir.Z());

        CP::THandle<CP::TClusterState> clusterState = cluster->GetState();
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                double v = clusterState->GetPositionCovariance(i,j);
                trackState->SetPositionCovariance(i,j,v);
                v = dirCov(i,j);
                trackState->SetDirectionCovariance(i,j,v);
            }
        }

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
    input->SetAlgorithmName("TClusterTrackFit");
    input->SetQuality(-logLikelihood);
    input->SetNDOF(trackDOF);

    return input;
}
