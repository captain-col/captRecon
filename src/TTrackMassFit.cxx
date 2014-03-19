#include "TTrackMassFit.hxx"
#include "TEnergyLoss.hxx"

#include <TReconTrack.hxx>
#include <TReconCluster.hxx>
#include <TReconNode.hxx>
#include <TRuntimeParameters.hxx>
#include <HEPUnits.hxx>

#include <Math/Factory.h>
#include <Math/Minimizer.h>
#include <Math/Functor.h>
#include <Math/IFunction.h>

#include <vector>

namespace CP {
    namespace MF {
        class TMassFitFunction;
    }
}

/// Return the likelihood for a particular set of mass fit parameters.
class CP::MF::TMassFitFunction : public ROOT::Math::IBaseFunctionMultiDim {
public:
    TMassFitFunction() : fForwardFit(true), 
                         fStoppingPenalty(1.0*unit::MeV) {
        fEnergyPerCharge = CP::TRuntimeParameters::Get().GetParameterD(
            "captRecon.energyPerCharge");
    }
    
    /// Calculate the function.
    double DoEval(const double *params) const;

    unsigned int NDim() const {return 3;}

    TMassFitFunction* Clone() const {return new TMassFitFunction(*this);}

    /// A flag to say which way the fit is going.  This is true if the fit
    /// goes from front to back (the stopping/exit point would be at the
    /// back).
    bool fForwardFit;

    /// A penalty to apply if the particle has a non-zero kinetic energy when
    /// it "stops".  This has units of MeV.  This is ignored if it is less
    /// than or equal to zero.
    double fStoppingPenalty;

    /// The amount of energy per measured ionization electron.  The charge
    /// conversion is set using the parameter file, but can be overridden once
    /// the object is created.
    double fEnergyPerCharge;

    /// The amount of deposited energy at this point.  This is in units of
    /// ionization electrons.
    std::vector<double> fEnergyDeposit;

    /// The uncertainty on each energy measurement.
    std::vector<double> fEnergySigma;

    /// The length over which the energy is deposited.
    std::vector<double> fLength;

};

double CP::MF::TMassFitFunction::DoEval(const double *params) const {
    double mass = params[0];
    double kinOffset = params[1];
    double energyScale = std::exp(params[2])*fEnergyPerCharge; 

    CP::TEnergyLoss energyLoss;

    // Set the range of energy deposits to fit over.  If there are enough
    // measurements, then exclude the first and last measurements since they
    // are oftened "confused" by overlapping tracks.
    int start = 0;
    int stop = fEnergyDeposit.size();
    if (fEnergyDeposit.size() > 5) {
        ++start;
        --stop;
    }

    double kinEnergy = kinOffset;
    double logLikelihood = 0.0;

    // Add a logarithmic penalty if the particle which should stop has a
    // offset energy.  If the particle is exiting (or otherwise, doesn't
    // stop), fStoppingPenalty should be negative.
    if (fStoppingPenalty > 0.0) {
        double v = kinOffset/fStoppingPenalty;
        logLikelihood += v;
    }

    double sFactor = 0.5;
    if (fForwardFit) {
        // Correct for not starting at the beginning...
        for (int i=0; i<start; ++i) {
            int index = fEnergyDeposit.size() - i - 1;
            double energyDeposit = energyScale*fEnergyDeposit[index];
            kinEnergy += energyDeposit;
        }
        for (int i=start; i<stop; ++i) {
            int index = fEnergyDeposit.size() - i - 1;
            double energyDeposit = energyScale*fEnergyDeposit[index];
            double energySigma = energyScale*fEnergySigma[index];
            double length = fLength[index];
            kinEnergy += sFactor*energyDeposit;
#ifdef PRINTIT
#undef PRINTIT
            double mostProbable = energyLoss.GetMostProbable(kinEnergy,
                                                             mass,
                                                             length);
            std::cout << "   " << index
                      << " " << kinEnergy
                      << " " << kinOffset
                      << " " << energyDeposit 
                      << " " << mostProbable 
                      << " " << energyDeposit/length
                      << " " << mostProbable/length
                      << " " << energySigma 
                      << std::endl;
#endif
            double likelihood 
                = energyLoss.GetDepositPDF(kinEnergy,
                                           mass,
                                           length,
                                           energyDeposit,
                                           energySigma);
            logLikelihood -= std::log(likelihood);
            kinEnergy += (1.0-sFactor)*energyDeposit;
        }
    }
    else {
        // Correct for not starting at the beginning...
        for (int i=0; i<start; ++i) {
            int index = i;
            double energyDeposit = energyScale*fEnergyDeposit[index];
            kinEnergy += energyDeposit;
        }
        for (int i=start; i<stop; ++i) {
            int index = i;
            double energyDeposit = energyScale*fEnergyDeposit[index];
            double energySigma = energyScale*fEnergySigma[index];
            double length = fLength[index];
            kinEnergy += sFactor*energyDeposit;
#ifdef PRINTIT
#undef PRINTIT
            double mostProbable = energyLoss.GetMostProbable(kinEnergy,
                                                             mass,
                                                             length);
            std::cout << "   " << index
                      << " " << kinEnergy
                      << " " << kinOffset
                      << " " << energyDeposit 
                      << " " << energySigma 
                      << " " << mostProbable 
                      << " " << energyDeposit/length
                      << std::endl;
#endif
            double likelihood 
                = energyLoss.GetDepositPDF(kinEnergy,
                                           mass,
                                           length,
                                           energyDeposit,
                                           energySigma);
            logLikelihood -= std::log(likelihood);
            kinEnergy += (1.0-sFactor)*energyDeposit;
        }
    }

    // Add a penalty if the offset gets too large.  This prevents particles
    // from climbing the Landau plateau.
    double mGam = energyLoss.GetMIPGamma();
    double r = kinOffset/mass;
    if (r > (mGam-1.0)) {
        double v = (kinOffset-(mGam-1.0)*mass);
        logLikelihood += 0.5*v*v;
    }

    // If the offset approaches MIP, the penalize the likelihood if the mass
    // is bigger that a muon.
    double penalty = (std::max(mass,105.0)-105.0)*r/(mGam-1);
    logLikelihood += penalty*penalty;

    penalty = params[2]/0.05;
    logLikelihood += penalty*penalty;

#ifdef PRINTIT
#undef PRINTIT
    std::cout << "logL: " << logLikelihood
              << " kinEnergy: " << kinEnergy/unit::MeV 
              << " mass: " << mass/unit::MeV
              << " offset: " << kinOffset/unit::MeV 
              << " scale: " << std::exp(params[2])
              << std::endl;
#endif

    return logLikelihood;
}


CP::TTrackMassFit::TTrackMassFit() {}
CP::TTrackMassFit::~TTrackMassFit() {}

namespace {
    void RunRootMinimizer(ROOT::Math::Minimizer& rootMinimizer,
                          double startMass, double startOffset,
                          double startScale) {

        rootMinimizer.SetLimitedVariable(0, "mass", startMass, 0.1,
                                         0.95*startMass, 1000.0*unit::MeV);

        rootMinimizer.SetLowerLimitedVariable(1, "kinOffset", startOffset, 0.1,
                                              -5.0*unit::MeV);
        
        rootMinimizer.SetVariable(2,"log(EScale)", std::log(startScale), 0.1);

        rootMinimizer.Minimize();
    }
}

CP::THandle<CP::TReconTrack> 
CP::TTrackMassFit::Apply(CP::THandle<CP::TReconTrack>& input) {
    TReconNodeContainer& nodes = input->GetNodes();
    
    if (nodes.size() < 4) {
        CaptError("Not enough nodes to fit.");
        return CP::THandle<CP::TReconTrack>();
    }

    /// Create the minimizer.
    CP::MF::TMassFitFunction fitFunction;
    std::auto_ptr<ROOT::Math::Minimizer> 
        rootMinimizer(ROOT::Math::Factory::CreateMinimizer("Minuit2"));
    rootMinimizer->SetFunction(fitFunction);
    rootMinimizer->SetErrorDef(1.0);

    // Get the energy deposition and cluster length as a function of position
    // along the track
    double totalEnergy = 0.0;
    for (TReconNodeContainer::iterator n = nodes.begin();
         n != nodes.end(); ++n) {
        CP::THandle<CP::TTrackState> state = (*n)->GetState();
        CP::THandle<CP::TReconCluster> object = (*n)->GetObject();
        TVector3 prevPos;
        if (n != nodes.begin()) {
            CP::THandle<CP::TTrackState> prevState = (*(n-1))->GetState();
            prevPos = 0.5*(prevState->GetPosition().Vect()
                           + state->GetPosition().Vect());
        }
        else {
            CP::THandle<CP::TTrackState> prevState = input->GetFront();
            prevPos = prevState->GetPosition().Vect();
        }
        TVector3 nextPos;
        if (n+1 != nodes.end()) {
            CP::THandle<CP::TTrackState> nextState = (*(n+1))->GetState();
            nextPos = 0.5*(nextState->GetPosition().Vect()
                           + state->GetPosition().Vect());
        }
        else {
            CP::THandle<CP::TTrackState> nextState = input->GetFront();
            nextPos = nextState->GetPosition().Vect();
        }
        double energy = state->GetEDeposit();
        double energySigma = std::sqrt(state->GetEDepositVariance());
        double length1 = (nextPos - prevPos).Mag();
        double length2 = 2.0*object->GetLongExtent();
        double length = std::min(length1,length2);
        totalEnergy += energy;
#ifdef PRINTIT
#undef PRINTIT
        std::cout << object->GetUniqueID()
                  << " " << 34.1*unit::eV*totalEnergy
                  << " " << 34.1*unit::eV*energy 
                  << " " << 34.1*unit::eV*energySigma
                  << " " << length
                  << " " << 34.1*unit::eV*unit::cm*energy/length << std::endl;
#endif
        fitFunction.fEnergyDeposit.push_back(energy);
        fitFunction.fEnergySigma.push_back(energySigma);
        fitFunction.fLength.push_back(length);
    }

    // Define the starting point of the fit.
    double startMass = 105.0*unit::MeV;
    double startOffset = 1.0*unit::MeV;
    double startScale = 1.0;

    fitFunction.fForwardFit = true;
    fitFunction.fStoppingPenalty = 3.0;
    RunRootMinimizer(*rootMinimizer, startMass,startOffset,startScale);

    double forwardMin = rootMinimizer->MinValue();
    // if (rootMinimizer->Status() != 0) forwardMin = 1E+200;
    double forwardMass = rootMinimizer->X()[0];
    double forwardMassSigma = rootMinimizer->Errors()[0];
    double forwardOffset = rootMinimizer->X()[1];
    double forwardScale = std::exp(rootMinimizer->X()[2]);
  
    CaptNamedInfo("TTrackMassFit", 
                  "forward UID: " << input->GetUniqueID()
                  << " size: " << nodes.size()
                  << " L: " << forwardMin
                  << " M: " << forwardMass
                  << " E: " << forwardOffset
                  << " S: " << forwardScale);

    double bestValue = forwardMin;
    int bestDirection = 1;
    double bestMass = forwardMass;
    double bestMassSigma = forwardMassSigma;
    double bestOffset = forwardOffset;
    double bestScale = forwardScale;

    fitFunction.fForwardFit = false;
    fitFunction.fStoppingPenalty = 3.0;
    RunRootMinimizer(*rootMinimizer, startMass,startOffset,startScale);

    double backwardMin = rootMinimizer->MinValue();
    // if (rootMinimizer->Status() != 0) backwardMin = 1E+200;
    double backwardMass = rootMinimizer->X()[0];
    double backwardMassSigma = rootMinimizer->Errors()[0];
    double backwardOffset = rootMinimizer->X()[1];
    double backwardScale = std::exp(rootMinimizer->X()[2]);
    CaptNamedInfo("TTrackMassFit", 
                  "backward UID: " << input->GetUniqueID()
                  << " size: " << nodes.size()
                  << " L: " << backwardMin
                  << " M: " << backwardMass
                  << " E: " << backwardOffset
                  << " S: " << backwardScale);

    // Find the best fit to the mass.
    if (bestValue > backwardMin) {
        bestValue = backwardMin;
        bestDirection = -1;
        bestMass = backwardMass;
        bestMassSigma = backwardMassSigma;
        bestOffset = backwardOffset;
        bestScale = backwardScale;
    }

    if (!(bestMassSigma > 0.01)) bestMassSigma = 0.01;

    CaptNamedLog("TTrackMassFit", 
                 "UID: " << input->GetUniqueID()
                 << " size: " << nodes.size()
                 << " L: " << bestValue
                 << " D: " << bestDirection
                 << " M: " << bestMass << " +/- " << bestMassSigma

                 << " E: " << bestOffset
                 << " S: " << bestScale);

    if (bestValue > 1E+100) return CP::THandle<CP::TReconTrack>();

    // Set the track node values.
    CP::THandle<CP::TTrackState> state = input->GetFront();
    state->SetMass(bestMass);
    state->SetMassVariance(bestMassSigma*bestMassSigma);

    state = input->GetBack();
    state->SetMass(bestMass);
    state->SetMassVariance(bestMassSigma*bestMassSigma);

    for (CP::TReconNodeContainer::iterator n = nodes.begin(); 
         n != nodes.end(); ++n) {
        state = (*n)->GetState();
        state->SetMass(bestMass);
        state->SetMassVariance(bestMassSigma*bestMassSigma);
    }

    if (bestDirection<0) input->ReverseTrack();

    return input;
}
