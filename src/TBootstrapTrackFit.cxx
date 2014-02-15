#include "TBootstrapTrackFit.hxx"

#include <TReconNode.hxx>
#include <TReconCluster.hxx>
#include <HEPUnits.hxx>
#include <HEPConstants.hxx>
#include <TCaptLog.hxx>

#include <ostreamTVector3.hxx>
#include <ostreamTLorentzVector.hxx>

#include <bfl/filter/bootstrapfilter.h>
#include <bfl/model/systemmodel.h>
#include <bfl/model/measurementmodel.h>
#include <bfl/pdf/conditionalpdf.h>
#include <bfl/pdf/gaussian.h>

#include <bfl/wrappers/rng/rng.h>

#include <TRandom.h>
#include <TDecompChol.h>

/// A "local" name space for the Bootstrap Track Fitter (BTF).  This is the
/// name space for all of the user classes needed to interface to the BFL.
/// These classes are only used in the TBootstrapTrackFit.cxx file.  This
/// could be placed in the anonymous namespace, but placing them in a named
/// one makes life a little easier.
namespace BTF {
    /// A class to update the prediction of the system state.
    class TSystemPDF;

    /// A class to calculate the probability of an expected measurement
    /// relative to the actual measurement.
    class TMeasurementPDF;

    /// A local wrapper around the BFL ColumnVector class.
    typedef MatrixWrapper::ColumnVector ColumnVector;

    /// An enumeration defining the state fields.
    enum StateDefinition {
        kXPos = 0, kYPos, kZPos,
        kXDir, kYDir, kZDir,
        kStateSize        
    };

    /// A function to fill the vector of prior states with samples from the
    /// prior PDF.  The prior PDF is based on a pair of nodes.  This generates
    /// the PDF assuming that the track is going from node1 to node2, and
    /// works for both forward and backward filtering.  However, for backward
    /// filtering, the prior is going to include a little bit of extra
    /// variance since it will add any multiple scattering that might have
    /// been introduced between node1 and node 2.  It's a small amount and as
    /// long as the track is long enough that the prior doesn't dominate, it
    /// is not a problem.  In anycase, if the track is short enough that the
    /// prior matters, then the fit and it's covariance is not reliable.
    void GeneratePrior(std::vector< BFL::Sample<BTF::ColumnVector> >& samples,
                       CP::THandle<CP::TReconNode> node1,
                       CP::THandle<CP::TReconNode> node2);

    /// This takes a covariance matrix that might have extreme correlations
    /// and "conditions" it by reducing the correlations until it is positive
    /// definite.  The corr parameter controls how much correlation is
    /// required.  This modifies the input matrix.
    void ConditionMatrix(TMatrixD& covar1, double corr = 0.01);
};

/// A class to update the prediction of the system state.
class BTF::TSystemPDF 
    : public BFL::ConditionalPdf<ColumnVector,ColumnVector> {
public:
    TSystemPDF(double momentumEstimate) 
        : BFL::ConditionalPdf<ColumnVector,ColumnVector>(kStateSize,1) {
        fMomentumEstimate = momentumEstimate;
    }
    virtual ~TSystemPDF() {}
    
    /// Set whether this is a forward or backward update.  Comment: The
    /// SampleFrom method is written so that it is independent of whether
    /// the class is used for forward or backward filtering.
    void SetBackwardUpdate(bool back=true) {fBackwardUpdate = back;}
    
    /// Set the measurement object (i.e. the cluster) associated with the
    /// next state.
    void SetMeasurement(CP::THandle<CP::TReconCluster> next) {
        fMeasurement = next;
    }
    
    /// The function to generate one updated state.  The previous state of
    /// the system is accessed through the ConditionalArgumentGet(0)
    /// method.  The updated state is set as the value of oneSample.  Note
    /// that the updated state is propagated to the expected state for the
    /// next measurement, and includes the effect of multiple scattering
    /// (and any other random proceses.  The final two arguments (method
    /// and args) are not used by the BFL, and always have the value of
    /// DEFAULT and NULL.
    virtual bool SampleFrom(BFL::Sample<ColumnVector>& oneSample,
                            int method = DEFAULT, 
                            void* args = NULL) const;
    
private:
    /// This is true if the PDF is being used for back-propagation.
    bool fBackwardUpdate;
    
    /// This is the measurement that will be associated with the state.
    CP::THandle<CP::TReconCluster> fMeasurement;

    /// The momentum estimate needed for multiple scattering.
    double fMomentumEstimate;
};

/// A class to calculate the probability of an expected measurement
/// relative to the actual measurement.
class BTF::TMeasurementPDF
    : public BFL::ConditionalPdf<ColumnVector,ColumnVector> {
public:
    TMeasurementPDF()
        : BFL::ConditionalPdf<ColumnVector,ColumnVector>(1,1) {}
    virtual ~TMeasurementPDF() {}
    
    /// Return the probability of a state based on a particular
    /// measurement.  The candidate state of the system is accessed
    /// through the ConditionalArgumentGet(0) method.  The return value is
    /// a specialized double that must be between 0.0 and 1.0 inclusive.
    virtual BFL::Probability
    ProbabilityGet(const ColumnVector& measurement) const;
    
    /// Get the probability for a position based on the current measurement.  
    virtual BFL::Probability
    ProbabilityGet(const TVector3& position)const;

    /// Set the measurement object (i.e. the cluster).  This also sets the
    /// value of the fGaussian field so that it describes the current
    /// measurement.
    void SetMeasurement(CP::THandle<CP::TReconCluster> next);
    
private:
    /// This is the measurement that will be associated with the state.
    CP::THandle<CP::TReconCluster> fMeasurement;

    /// The position of the measurement
    TVector3 fPosition;

    /// The error matrix (inverted covariance) of the measurement.
    TMatrixD fErrorMat;

    /// The gaussian prefactor.
    double fGaussianConstant;
};

CP::TBootstrapTrackFit::TBootstrapTrackFit(int trials) : fTrials(trials) { }
CP::TBootstrapTrackFit::~TBootstrapTrackFit() {}

double
CP::TBootstrapTrackFit::EstimateRange(const CP::TReconNodeContainer& nodes) { 
    CP::THandle<CP::TReconCluster> lastCluster;
    double range = 0.0;
    for (CP::TReconNodeContainer::const_iterator n = nodes.begin();
         n != nodes.end(); ++n) { 
        CP::THandle<CP::TReconCluster> cluster = (*n)->GetObject();
        if (!cluster) {
            throw EBootstrapMissingNodeObject();
        }
        if (lastCluster) {
            range += (lastCluster->GetPosition().Vect() 
                      - cluster->GetPosition().Vect()).Mag();
        }
        lastCluster = cluster;
    }
    return range;
}

double
CP::TBootstrapTrackFit::EstimateMomentum(const CP::TReconNodeContainer& nodes) {
    double range = EstimateRange(nodes);
    double dEdX = 0.2*unit::MeV/unit::mm;
    double offset = 40*unit::MeV;
    double mass = 105*unit::MeV;
    double energy = range*dEdX + offset + mass;
    double momentum = std::sqrt(energy*energy - mass*mass);
    return momentum;
}

CP::THandle<CP::TReconTrack>
CP::TBootstrapTrackFit::Apply(CP::THandle<CP::TReconTrack>& input) {

    TReconNodeContainer& nodes = input->GetNodes();
    if (nodes.size() < 2) {
        CaptError("Not enough nodes to fit.");
        return CP::THandle<CP::TReconTrack>();
    }

    double range = EstimateRange(nodes);
    double momentum = EstimateMomentum(nodes);

    CaptLog("Bootstrap " << input->GetUniqueID() 
            << " w/ " << nodes.size() << " nodes"
            << "    range: " << range/unit::mm << " mm"
            << "    rough momentum: " << momentum << " MeV/c");

    // Define the number of states to use in the particle filter.
    const int numSamples = fTrials;

    // Define the size of the state.
    const int stateSize = BTF::kStateSize;

    BTF::TSystemPDF systemPDF(momentum);
    BFL::SystemModel<BTF::ColumnVector> systemModel(&systemPDF);
    BTF::TMeasurementPDF measurementPDF;
    BFL::MeasurementModel<BTF::ColumnVector, BTF::ColumnVector> 
        measurementModel(&measurementPDF);
    BTF::ColumnVector measurement(1);

    /////////////////////////////////////////////////////////////////////
    // Forward filtering.
    /////////////////////////////////////////////////////////////////////
    systemPDF.SetBackwardUpdate(false);

    /////////////////////////////////////////////////////////////////////
    // Generate the prior samples that are used to start the filter.  The
    // prior samples are based on the first two nodes and are
    // positioned just "upstream" of the first node.
    /////////////////////////////////////////////////////////////////////
    CP::THandle<CP::TReconNode> firstNode = nodes[0];
    CP::THandle<CP::TReconNode> otherNode = nodes[1];
    for (std::size_t i=1; i<nodes.size()/3; ++i) {
        CP::THandle<CP::TReconCluster> cluster1 = firstNode->GetObject();
        CP::THandle<CP::TReconCluster> cluster2 = nodes[i]->GetObject();
        double diff = (cluster1->GetPosition().Vect() 
                       - cluster2->GetPosition().Vect()).Mag();
        if (diff > 2*unit::cm) break;
        otherNode = nodes[i];
    }
    std::vector< BFL::Sample<BTF::ColumnVector> > priorSamples(numSamples);
    BTF::GeneratePrior(priorSamples, firstNode, otherNode);
    BFL::MCPdf<BTF::ColumnVector> forwardPrior(numSamples,stateSize);
    forwardPrior.ListOfSamplesSet(priorSamples);

    BFL::BootstrapFilter<BTF::ColumnVector, BTF::ColumnVector> 
        forwardFilter(&forwardPrior, 0, numSamples/4.0);

    /////////////////////////////////////////////////////////////////////
    /// Estimate the state at each node.
    /////////////////////////////////////////////////////////////////////
    for (std::size_t step = 0; step < nodes.size(); ++step) {
        measurement[0] = step;
        CP::THandle<CP::TReconNode> node = nodes[step];

        // Get the measurement from the node.
        CP::THandle<CP::TTrackState> trackState = node->GetState();
        CP::THandle<CP::TReconCluster> cluster = node->GetObject();
        systemPDF.SetMeasurement(cluster);
        measurementPDF.SetMeasurement(cluster);

        // Update the filter
        forwardFilter.Update(&systemModel, &measurementModel, measurement);
        
        // Get the updated state and it's covariance.  
        BFL::MCPdf<BTF::ColumnVector>* posterior = forwardFilter.PostGet(); 

        // Use the current expectation value and covariance to give temporary
        // values to the track nodes.
        BTF::ColumnVector ev = posterior->ExpectedValueGet();
        MatrixWrapper::SymmetricMatrix cv = posterior->CovarianceGet();

        // Set the values and covariances.
        trackState->SetEDeposit(cluster->GetEDeposit());
        trackState->SetPosition(ev[BTF::kXPos],ev[BTF::kYPos],ev[BTF::kZPos],
                                cluster->GetPosition().T());
        trackState->SetDirection(ev[BTF::kXDir],ev[BTF::kYDir],ev[BTF::kZDir]);

        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                double v = cv(BTF::kXPos+i+1, BTF::kXPos+j+1);
                trackState->SetPositionCovariance(i,j,v);
                v = cv(BTF::kXDir+i+1, BTF::kXDir+j+1);
                trackState->SetDirectionCovariance(i,j,v);
            }
        }

        // This is the first step, so give it a large variance.  Otherwise the
        // prior has too much weight.
        if (step < 1) {
            for (int i=0; i<3; ++i) {
                trackState->SetPositionCovariance(i,i,10*unit::cm);
                trackState->SetDirectionCovariance(i,i,1.0);
            }
        }

        CaptNamedVerbose("BFL","Forward " 
                         << trackState->GetPosition().Vect()
                         << " -> " << trackState->GetDirection()
                         << " diff " 
                         << (trackState->GetPosition().Vect()
                             - cluster->GetPosition().Vect()).Mag());
        
    }

    /////////////////////////////////////////////////////////////////////
    // Backward filtering.
    /////////////////////////////////////////////////////////////////////
    systemPDF.SetBackwardUpdate(true);

    /////////////////////////////////////////////////////////////////////
    // Generate the prior samples that are used to start the filter.  The
    // prior samples are based on the first two nodes and are
    // positioned just "upstream" of the first node.
    /////////////////////////////////////////////////////////////////////
    firstNode = nodes[nodes.size()-1];
    otherNode = nodes[nodes.size()-2];
    for (std::size_t i=2; i<nodes.size()/3; ++i) {
        CP::THandle<CP::TReconCluster> cluster1 = firstNode->GetObject();
        CP::THandle<CP::TReconCluster> cluster2 
            = nodes[nodes.size()-i]->GetObject();
        double diff = (cluster1->GetPosition().Vect() 
                       - cluster2->GetPosition().Vect()).Mag();
        if (diff > 2*unit::cm) break;
        otherNode = nodes[nodes.size()-i];
    }
    BTF::GeneratePrior(priorSamples, otherNode, firstNode);
    BFL::MCPdf<BTF::ColumnVector> backwardPrior(numSamples,stateSize);
    backwardPrior.ListOfSamplesSet(priorSamples);

    BFL::BootstrapFilter<BTF::ColumnVector, BTF::ColumnVector> 
        backwardFilter(&backwardPrior, 0, numSamples/4.0);

    // Calculate the goodness of fit during the backfilter.
    double logLikelihood = 0.0;

    /////////////////////////////////////////////////////////////////////
    /// Estimate the state at each node.
    /////////////////////////////////////////////////////////////////////
    for (int step = (int) nodes.size()-1; 0 <= step; --step) {
        CaptNamedDebug("BFL","Start backward step " << step);
        measurement[0] = step;
        CP::THandle<CP::TReconNode> node = nodes[step];
        if (!node) {
            CaptError("Missing node in step " << step);
            continue;
        }

        // Get the measurement from the node.
        CP::THandle<CP::TTrackState> trackState = node->GetState();
        if (!trackState) {
            CaptError("Missing track state in step " << step);
            continue;
        }

        CP::THandle<CP::TReconCluster> cluster = node->GetObject();
        if (!cluster) {
            CaptError("Missing cluster in step " << step);
            continue;
        }

        systemPDF.SetMeasurement(cluster);
        measurementPDF.SetMeasurement(cluster);

        // Update the filter
        backwardFilter.Update(&systemModel, &measurementModel, measurement);
        
        // Get the updated state and it's covariance.
        BFL::MCPdf<BTF::ColumnVector>* posterior = backwardFilter.PostGet();
        BTF::ColumnVector ev = posterior->ExpectedValueGet();
        MatrixWrapper::SymmetricMatrix cv = posterior->CovarianceGet();

        // Set the values and covariances based on the back filtering, but
        // drop the correlations between position and direction (it makes life
        // simpler, and makes the calculation more stable).  

        // First check that the covariances are OK.
        double oldCov = 0.0;
        double newCov = 0.0;
        for (int i=0; i<3; ++i) {
            oldCov += trackState->GetPositionCovariance(i,i);
            newCov += cv(BTF::kXPos+i+1,BTF::kXPos+i+1);
        }

        if (oldCov < 1E-6 || newCov < 1E-6) {
            CaptNamedError("BFL", "Fit failed");
            return CP::THandle<CP::TReconTrack>();
        }

        // Save the forward filtering position, direction and variances.
        TVector3 posF(trackState->GetPosition().Vect());
        TVector3 vPosF(trackState->GetPositionVariance().Vect());
        TVector3 dirF(trackState->GetDirection());
        double vDirF = 0.0;
        for (int i=0; i<3; ++i) {
            vDirF += trackState->GetDirectionCovariance(i,i);
            // Turn the variance into a sigma, and make sure that zeros have
            // no weight.
            vPosF[i] = std::sqrt(vPosF[i]);
            if (vPosF[i] < 0.001*unit::mm) vPosF[i] = 1000*unit::meter;
        }

        // Save the backward filtering position, direction and variances
        TVector3 posB(ev[BTF::kXPos],ev[BTF::kYPos],ev[BTF::kZPos]);
        TVector3 vPosB(cv(BTF::kXPos+1, BTF::kXPos+1),
                       cv(BTF::kYPos+1, BTF::kYPos+1),
                       cv(BTF::kZPos+1, BTF::kZPos+1));
        TVector3 dirB(ev[BTF::kXDir],ev[BTF::kYDir],ev[BTF::kZDir]);
        double vDirB = 0.0;
        for (int i=0; i<3; ++i) {
            vDirB += cv(BTF::kXDir+i+1, BTF::kXDir+i+1);
            // Turn the variance into a sigma, and make sure that zeros have
            // no weight.
            vPosB[i] = std::sqrt(vPosB[i]);
            if (vPosB[i] < 0.001*unit::mm) vPosB[i] = 1000*unit::meter;
        }

        // Find the average position.
        double xx = posF.X()/vPosF.X() + posB.X()/vPosB.X();
        xx = xx/(1.0/vPosF.X() + 1.0/vPosB.X());
        double yy = posF.Y()/vPosF.Y() + posB.Y()/vPosB.Y();
        yy = yy/(1.0/vPosF.Y() + 1.0/vPosB.Y());
        double zz = posF.Z()/vPosF.Z() + posB.Z()/vPosB.Z();
        zz = zz/(1.0/vPosF.Z() + 1.0/vPosB.Z());
        
        // Find the average direction.
        TVector3 dir = (1.0/vDirF)*dirF + (1.0/vDirB)*dirB;
        dir = dir.Unit();

        // Set the final state.
        trackState->SetPosition(xx, yy, zz, cluster->GetPosition().T());
        trackState->SetDirection(dir.X(), dir.Y(), dir.Z());

        TMatrixD pcov1(3,3);
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                pcov1(i,j) = trackState->GetPositionCovariance(i,j);
            }
        }
        pcov1.InvertFast();

        TMatrixD dcov1(3,3);
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                dcov1(i,j) = trackState->GetDirectionCovariance(i,j);
            }
        }
        dcov1.InvertFast();

        TMatrixD pcov2(3,3);
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                pcov2(i,j) = cv(BTF::kXPos+i+1, BTF::kXPos+j+1);
                if (step == (int) nodes.size()-1 && i == j) {
                    pcov2(i,j) = 10*unit::cm;
                }
            }
        }
        pcov2.InvertFast();

        TMatrixD dcov2(3,3);
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                dcov2(i,j) = cv(BTF::kXDir+i+1, BTF::kXDir+j+1);
                if (step == (int) nodes.size()-1 && i == j) dcov2(i,j) = 1.0;
            }
        }
        dcov2.InvertFast();
        
        TMatrixD pcov(3,3);
        pcov = pcov1 + pcov2;
        pcov.InvertFast();
        TMatrixD dcov(3,3);
        dcov = dcov1 + dcov2;
        dcov.InvertFast();

        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                trackState->SetPositionCovariance(i,j,pcov(i,j));
                trackState->SetDirectionCovariance(i,j,dcov(i,j));
            }
        }

        CaptNamedVerbose("BFL","Backward " 
                         << trackState->GetPosition().Vect()
                         << " -> " << trackState->GetDirection()
                         << " diff " 
                         << (trackState->GetPosition().Vect()
                             - cluster->GetPosition().Vect()).Mag());
        
        double r = measurementPDF.ProbabilityGet(
            trackState->GetPosition().Vect());
        if (r>0) logLikelihood += std::log(r);
    }

    double energyDeposit = 0.0;
    for (TReconNodeContainer::iterator n = nodes.begin();
         n != nodes.end(); ++n) {
        CP::THandle<CP::TTrackState> state = (*n)->GetState();
        energyDeposit += state->GetEDeposit();
    }

    // Set the front state of the track.
    CP::THandle<CP::TTrackState> frontState = input->GetFront();
    CP::THandle<CP::TTrackState> firstNodeState = nodes.front()->GetState();
    *frontState = *firstNodeState;
    frontState->SetEDeposit(energyDeposit);

    if (!nodes.front()->GetObject()) {
        // Make sure there is an object, and there had *better* be one or
        // something is going horribly wrong.
        CaptError("Missing object for the front node");
        throw EBootstrapMissingNodeObject();
    }

    // Now move the front state upstream to the position of the first hit.
    // Notice that the front state is not at the same location as the first
    // node.  
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

    if (!nodes.back()->GetObject()) {
        // Make sure there is an object, and there had *better* be one or
        // something is going horribly wrong.
        CaptError("Missing object for the front node");
        throw EBootstrapMissingNodeObject();
    }

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

    int trackDOF = 3*nodes.size() - 6;
    input->SetStatus(TReconBase::kSuccess);
    input->SetAlgorithmName("TBootstrapTrackFit");
    input->SetQuality(-logLikelihood);
    input->SetNDOF(trackDOF);

    if (backDir * frontDir < 0.0) {
        CaptError("Front and back in opposite directions"
                  << " front " << frontDir
                  << " back " << backDir);
    }

    CaptLog("Bootstrap filter finished.  Quality: " 
            << - logLikelihood 
            << " / " << trackDOF);

    return input;
}

bool BTF::TSystemPDF::SampleFrom(BFL::Sample<ColumnVector>& oneSample,
                                 int method, void* args) const {
    // Get the previous state.
    ColumnVector state = ConditionalArgumentGet(0);

    TVector3 pos(state[BTF::kXPos], state[BTF::kYPos], state[BTF::kZPos]);
    TVector3 dir(state[BTF::kXDir], state[BTF::kYDir], state[BTF::kZDir]);

    for (int i=0; i<3; ++i) {
        if (!std::isfinite(pos[i])) {
            CaptError("Error in position " << pos);
            break;
        }
    }
            
    for (int i=0; i<3; ++i) {
        if (!std::isfinite(dir[i])) {
            CaptError("Error in direction " << dir);
            break;
        }
    }

    if (dir.Mag() < 0.1) {
            CaptError("Direction is short " << dir);
    }

    // Find the distance to the closest approach to the next measurement.
    TVector3 diff1 = fMeasurement->GetPosition().Vect() - pos;
    double dist1 = diff1*dir;

    // Find the amount of scattering. This can be found based on the
    // theoretical multiple scattering, or based on the deviation of a set of
    // nodes from a straight line.  The amount of scattering isn't super
    // crucial, but does determine how many nodes contribute to the state
    // information since the scatter will erase the effect of very distant
    // nodes.
    double posScatter = 0.0;
    double dirScatter = 0.0;
#define THEORY_SCATTER
#ifdef THEORY_SCATTER
    // Note that this is only valid for muon like tracks.  The values are
    // taken from the PDG.
    double radLen = 14*unit::cm; // For liquid argon.
    double X = std::abs(dist1)/radLen;
    double P = fMomentumEstimate/2.0; // HACK!  Assume have total momentum
    if (X < 0.001) X = 0.001;
    // Set the minimum amount of scattering.  This isn't very physical, but
    // there needs to be scattering or the fit doesn't work.
    dirScatter = (1.0+0.038*std::log(X))*std::sqrt(X)*(13.6*unit::MeV)/(P);
    posScatter = dist1*dirScatter;
#endif

    // Generate a random amount of multiple scattering and use it to update
    // the position and direction.  Make sure that the direction is
    // normalized.  I'm approximating actual fluctuations, so this will
    // slightly overestimate the position variance.  However, since this is
    // "track" fitting and not "particle" fitting, we don't have a PID
    // hypothesis yet, so the scattering isn't right anyway.

    if (dirScatter > 0.0) {
        double scatter = gRandom->Gaus(0.0, dirScatter);
        TVector3 ortho = scatter * (dir.Orthogonal().Unit());
        double phi = gRandom->Uniform(unit::twopi);
        ortho.Rotate(phi,dir);
        dir = (dir + ortho).Unit();
        for (int i=0; i<3; ++i) {
            pos[i] = pos[i] + gRandom->Gaus(0.0, posScatter);
        }
    }

    // Find the distance to the closest approach to the next measurement using
    // the updated direction.
    TVector3 diff2 = fMeasurement->GetPosition().Vect() - pos;
    double dist2 = diff2*dir;
    pos = pos + dist2*dir;

    for (int i=0; i<3; ++i) {
        if (!std::isfinite(pos[i])) {
            CaptError("Error in position " << pos);
            break;
        }
    }
            
    for (int i=0; i<3; ++i) {
        if (!std::isfinite(dir[i])) {
            CaptError("Error in direction " << dir);
            break;
        }
    }

    if (dir.Mag() < 0.1) {
            CaptError("Direction is short " << dir);
    }

    // Update the state.
    state[BTF::kXPos] = pos.X();
    state[BTF::kYPos] = pos.Y();
    state[BTF::kZPos] = pos.Z();
    state[BTF::kXDir] = dir.X();
    state[BTF::kYDir] = dir.Y();
    state[BTF::kZDir] = dir.Z();

    // Return the results in oneSample
    oneSample.ValueSet(state);

    return true;
}

BFL::Probability
BTF::TMeasurementPDF::ProbabilityGet(const ColumnVector& dummy) const {
    ColumnVector state = ConditionalArgumentGet(0);
    TVector3  expected(state[BTF::kXPos],state[BTF::kYPos],state[BTF::kZPos]);
    return ProbabilityGet(expected);
}

BFL::Probability
BTF::TMeasurementPDF::ProbabilityGet(const TVector3& expected) const {
    TVector3 difference = expected - fPosition;
    double X2 = 0.0;
    for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
            X2 += difference(i)*fErrorMat(i,j)*difference(j);
        }
    }
    double v = std::exp(-0.5*X2)*fGaussianConstant;
    if (!std::isfinite(v)) {
        CaptError("Invalid probability " << expected 
                  << " " << fPosition
                  << " " << X2 
                  << " " << fGaussianConstant
                  << " " << v);
    }
    return v;
}

void BTF::TMeasurementPDF::SetMeasurement(CP::THandle<CP::TReconCluster> next) {
    fMeasurement = next;
    CP::THandle<CP::TClusterState> state = fMeasurement->GetState();

    fPosition = state->GetPosition().Vect();
    fErrorMat.ResizeTo(3,3);
    for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
            double v = state->GetPositionCovariance(i,j);
            // Fix spurious accuracy for the Z axis.  The extra accuracy
            // is introduced by just taking the average of all of the
            // contributing hit positions.  However, the Z positions of
            // all of the hits are the same when the clusters are
            // constructed.  The intrinsic Z accuracy is about 0.5 mm.
            if (i == j) v = std::max(0.25*unit::mm*unit::mm, v);
            fErrorMat(i,j) = v;
        }
    }
    BTF::ConditionMatrix(fErrorMat);

    // The matrix has been filled with the covariance, not turn it into the
    // error matrix.
    double det;
    fErrorMat.InvertFast(&det);

    det = det*unit::pi2*unit::pi2*unit::pi2;
    fGaussianConstant = 1.0/sqrt(det);
}

void BTF::ConditionMatrix(TMatrixD& covar, double correlations) {
    TMatrixD mat(3,3);

    double scale = 1.0;
    for (int trials = 0; trials < 20; ++trials) {
        double var = 0.0;
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                double v = covar(i,j);
                if (i == j) var += v*v;
                // Iteratively reduce the amount of correlation between the
                // axes until the error matrix has a positive determinant.
                if (i != j) v *= scale;
                mat(i,j) = v;
            }
        }
        // The matrix has been filled with the covariance, not turn it into the
        // error matrix.
        double det;
        mat.InvertFast(&det);
        double corr = det/var;
        if (corr > correlations) break;
        scale = scale * 0.3;
    }
    
    for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
            if (i==j) continue;
            covar(i,j) *= scale;
        }
    }
}

void BTF::GeneratePrior(std::vector< BFL::Sample<BTF::ColumnVector> >& samples,
                        CP::THandle<CP::TReconNode> node1,
                        CP::THandle<CP::TReconNode> node2) {
    CP::THandle<CP::TReconCluster> cluster1 = node1->GetObject();
    CP::THandle<CP::TReconCluster> cluster2 = node2->GetObject();
    if (!cluster1 || !cluster2) {
        CaptError("Invalid objects (not TReconCluster) for BTF::GeneratePrior");
        std::exit(1);
    }

    CP::THandle<CP::TClusterState> state1 = cluster1->GetState();
    CP::THandle<CP::TClusterState> state2 = cluster2->GetState();
    if (!state1 || !state2) {
        CaptError("Invalid states (not TClusterState) for BTF::GeneratePrior");
        std::exit(1);
    }

    // Get the positions from both nodes.
    TVector3 pos1 = state1->GetPosition().Vect();
    TVector3 pos2 = state2->GetPosition().Vect();

    // Get the position covariances from both nodes.
    TMatrixD mat1(3,3);
    TMatrixD mat2(3,3);
    for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
            mat1(i,j) = state1->GetPositionCovariance(i,j);
            mat2(i,j) = state2->GetPositionCovariance(i,j);
        }
    }

    // Do the cholesky decomposition of the covariance matricies
    BTF::ConditionMatrix(mat1);
    TDecompChol decomp1(mat1);

    BTF::ConditionMatrix(mat2);
    TDecompChol decomp2(mat2);
    decomp1.Decompose();
    decomp2.Decompose();
    
    // And generate states with the right correllations. 
    BTF::ColumnVector sample(BTF::kStateSize);
    for (std::vector< BFL::Sample<BTF::ColumnVector> >::iterator 
             s = samples.begin();
         s != samples.end(); ++s) {
        TVector3 rdm1(gRandom->Gaus(), gRandom->Gaus(), gRandom->Gaus());
        TVector3 rdm2(gRandom->Gaus(), gRandom->Gaus(), gRandom->Gaus());
        rdm1 = decomp1.GetU()*rdm1 + pos1;
        rdm2 = decomp2.GetU()*rdm2 + pos2;
        TVector3 dir = (rdm2-rdm1).Unit();
        sample[BTF::kXPos] = rdm1.X();
        sample[BTF::kYPos] = rdm1.Y();
        sample[BTF::kZPos] = rdm1.Z();
        sample[BTF::kXDir] = dir.X();
        sample[BTF::kYDir] = dir.Y();
        sample[BTF::kZDir] = dir.Z();
        s->ValueSet(sample);
    }

};
