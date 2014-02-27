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

    unsigned int NDim() const {return 2;}

    TMassFitFunction* Clone() const {return new TMassFitFunction(*this);}

    /// A flag to say which way the fit is going.  This is true if the fit
    /// goes from front to back (the stopping/exit point would be at the
    /// back).
    bool fForwardFit;

    /// A penalty to apply if the particle has a non-zero kinetic energy when
    /// it "stops".  This has units of MeV.  This is ignored if it is less
    /// than or equal to zero.
    double fStoppingPenalty;

    /// The amount of energy per measured ionization electron.
    double fEnergyPerCharge;

    /// The amount of deposited energy at this point.  This is set using the
    /// parameter file, but can be overridden once the object is created.
    std::vector<double> fEnergyDeposit;

    /// The length over which the energy is deposited.
    std::vector<double> fLength;

};

double CP::MF::TMassFitFunction::DoEval(const double *params) const {
    double mass = std::exp(params[0]);
    double kinOffset = params[1];
    double energyScale = fEnergyPerCharge; 

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
    if (fStoppingPenalty > 0.0) {
        double v = kinOffset/fStoppingPenalty;
        logLikelihood += 0.5*v*v;
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
            double length = fLength[index];
            kinEnergy += sFactor*energyDeposit;
#ifdef PRINTIT
#undef PRINTIT
            double mostProbable = energyLoss.GetMostProbable(kinEnergy,
                                                             mass,
                                                             length);
            std::cout << "   " << index
                      << " " << kinEnergy
                      << " " << energyDeposit 
                      << " " << mostProbable 
                      << " " << energyDeposit/length
                      << std::endl;
#endif
            double likelihood 
                = energyLoss.GetDepositPDF(energyDeposit,
                                           kinEnergy,
                                           mass,
                                           length);
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
            double length = fLength[index];
            kinEnergy += sFactor*energyDeposit;
#ifdef PRINTIT
#undef PRINTIT
            double mostProbable = energyLoss.GetMostProbable(kinEnergy,
                                                             mass,
                                                             length);
            std::cout << "   " << index
                      << " " << kinEnergy
                      << " " << energyDeposit 
                      << " " << mostProbable 
                      << " " << energyDeposit/length
                      << std::endl;
#endif
            double likelihood 
                = energyLoss.GetDepositPDF(energyDeposit, 
                                           kinEnergy,
                                           mass,
                                           length);
            logLikelihood -= std::log(likelihood);
            kinEnergy += (1.0-sFactor)*energyDeposit;
        }
    }

#ifdef PRINTIT
#undef PRINTIT
    std::cout << "kinEnergy: " << kinEnergy/unit::MeV 
              << "  mass: " << mass
              << "  logL: " << logLikelihood << std::endl;
#endif

    return logLikelihood;
}


CP::TTrackMassFit::TTrackMassFit() {}
CP::TTrackMassFit::~TTrackMassFit() {}

CP::THandle<CP::TReconTrack> 
CP::TTrackMassFit::Apply(CP::THandle<CP::TReconTrack>& input) {
    TReconNodeContainer& nodes = input->GetNodes();
    
    if (nodes.size() < 4) {
        CaptError("Not enough nodes to fit.");
        return input;
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
        double length1 = (nextPos - prevPos).Mag();
        double length2 = 2.0*object->GetLongExtent();
        double length = std::min(length1,length2);
        totalEnergy += energy;
#ifdef PRINTIT
#undef PRINTIT
        std::cout << object->GetUniqueID()
                  << " " << 34.1*unit::eV*totalEnergy
                  << " " << 34.1*unit::eV*energy 
                  << " " << length
                  << " " << 34.1*unit::eV*unit::cm*energy/length << std::endl;
#endif
        fitFunction.fEnergyDeposit.push_back(energy);
        fitFunction.fLength.push_back(length);
    }

    fitFunction.fForwardFit = true;
    rootMinimizer->SetVariable(0,"log(Mass)", std::log(105.0*unit::MeV), 0.1);
    rootMinimizer->SetVariable(1,"kinOffset", 0.0*unit::MeV, 10*unit::MeV);
    rootMinimizer->Minimize();
    double forwardMin = rootMinimizer->MinValue();
    double forwardMass = std::exp(rootMinimizer->X()[0]);
    double forwardMassSigma = rootMinimizer->Errors()[0];
    double forwardOffset = std::abs(rootMinimizer->X()[1]);

    CaptNamedInfo("TTrackMassFit", 
                 "forward UID: " << input->GetUniqueID()
                 << " size: " << nodes.size()
                 << " L: " << forwardMin
                 << " M: " << forwardMass
                 << " E: " << forwardOffset);

    fitFunction.fForwardFit = false;
    rootMinimizer->SetVariable(0,"log(Mass)", std::log(105.0*unit::MeV), 0.1);
    rootMinimizer->SetVariable(1,"kinOffset", 0.0*unit::MeV, 10*unit::MeV);
    rootMinimizer->Minimize();
    double backwardMin = rootMinimizer->MinValue();
    double backwardMass = std::exp(rootMinimizer->X()[0]);
    double backwardMassSigma = rootMinimizer->Errors()[0];
    double backwardOffset = std::abs(rootMinimizer->X()[1]);

    CaptNamedInfo("TTrackMassFit", 
                 "backward UID: " << input->GetUniqueID()
                 << " size: " << nodes.size()
                 << " L: " << backwardMin
                 << " M: " << backwardMass
                 << " E: " << backwardOffset);

    double bestValue = forwardMin;
    int bestDirection = 1;
    double bestMass = forwardMass;
    double bestMassSigma = forwardMassSigma;
    double bestOffset = forwardOffset;
    
    // Find the best fit to the mass.
    if (bestValue > backwardMin) {
        bestValue = backwardMin;
        bestDirection = -1;
        bestMass = backwardMass;
        bestMassSigma = backwardMassSigma;
        bestOffset = backwardOffset;
    }

    CaptNamedLog("TTrackMassFit", 
                 "UID: " << input->GetUniqueID()
                 << " size: " << nodes.size()
                 << " L: " << bestValue
                 << " D: " << bestDirection
                 << " M: " << bestMass
                 << " S: " << bestMassSigma
                 << " E: " << bestOffset);

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
