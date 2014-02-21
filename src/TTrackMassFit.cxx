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
                         fStoppingPenalty(1.0*unit::MeV),
                         fScaleConstraint(0.01) {
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

    /// The amount of energy per measured ionization electron.
    double fEnergyPerCharge;

    /// A constraint on the scale parameter.  The scale give the fit a little
    /// bit of wiggle room round the nominal energy scale, and this paramter
    /// keeps the value close to 1.0.
    double fScaleConstraint;

    /// The amount of deposited energy at this point.  This is set using the
    /// parameter file, but can be overridden once the object is created.
    std::vector<double> fEnergyDeposit;

    /// The length over which the energy is deposited.
    std::vector<double> fLength;

};

double CP::MF::TMassFitFunction::DoEval(const double *params) const {
    double mass = std::exp(params[0]);
    double scale = std::exp(params[1]);
    double kinOffset = params[2];
    double energyScale = scale*fEnergyPerCharge; 

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

    if (fScaleConstraint > 1E-6) {
        double v = (scale-1.0)/fScaleConstraint;
        logLikelihood += 0.5*v*v;
    }

    if (fForwardFit) {
        for (int i=start; i<stop; ++i) {
            int index = fEnergyDeposit.size() - i - 1;
            double energyDeposit = energyScale*fEnergyDeposit[index];
            double length = fLength[index];
            kinEnergy += energyDeposit;
#define PRINTIT
#ifdef PRINTIT
            double mostProbable = energyLoss.GetMostProbable(kinEnergy,
                                                             mass,
                                                             length);
            std::cout << "   " << kinEnergy
                      << " " << energyDeposit 
                      << " " << mostProbable 
                      << " " << energyDeposit/length
                      << std::endl;
#endif
            double likelihood 
                = energyLoss.GetDepositPDF(energyDeposit, kinEnergy,
                                           mass, length);
            logLikelihood -= std::log(likelihood);
        }
    }
    else {
        for (int i=start; i<stop; ++i) {
            int index = i;
            double energyDeposit = energyScale*fEnergyDeposit[index];
            double length = fLength[index];
            kinEnergy += energyDeposit;
#ifdef PRINTIT
            double mostProbable = energyLoss.GetMostProbable(kinEnergy,
                                                             mass,
                                                             length);
            std::cout << "   " << kinEnergy
                      << " " << energyDeposit 
                      << " " << mostProbable 
                      << " " << energyDeposit/length
                      << std::endl;
#endif
            double likelihood 
                = energyLoss.GetDepositPDF(energyDeposit, kinEnergy,
                                           mass, length);
            logLikelihood -= std::log(likelihood);
        }
    }

#ifdef PRINTIT
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
    
    if (nodes.size() < 2) {
        CaptError("Not enough nodes to fit.");
        return CP::THandle<CP::TReconTrack>();
    }

    /// Create the minimizer.
    CP::MF::TMassFitFunction fitFunction;
    std::auto_ptr<ROOT::Math::Minimizer> 
        rootMinimizer(ROOT::Math::Factory::CreateMinimizer("Minuit2"));
    rootMinimizer->SetFunction(fitFunction);
    rootMinimizer->SetErrorDef(0.5);

    // Get the energy deposition and cluster length as a function of position
    // along the track
    for (TReconNodeContainer::iterator n = nodes.begin();
         n != nodes.end(); ++n) {
        CP::THandle<CP::TTrackState> state = (*n)->GetState();
        CP::THandle<CP::TReconCluster> object = (*n)->GetObject();
        double energy = state->GetEDeposit();
        double length = 2.0*object->GetLongExtent();
        fitFunction.fEnergyDeposit.push_back(energy);
        fitFunction.fLength.push_back(length);
    }

    fitFunction.fForwardFit = true;
    rootMinimizer->SetVariable(0,"log(Mass)", std::log(105.0*unit::MeV), 0.1);
    rootMinimizer->SetVariable(1,"log(Scale)", 0.0, 0.1);
    rootMinimizer->SetVariable(2,"kinOffset", 0.0*unit::MeV, 10*unit::MeV);
    rootMinimizer->Minimize();
    double forwardMin = rootMinimizer->MinValue();
    double forwardMass = std::exp(rootMinimizer->X()[0]);
    double forwardScale = std::exp(rootMinimizer->X()[1]);
    double forwardOffset = std::abs(rootMinimizer->X()[2]);

    CaptNamedLog("TTrackMassFit", 
                 "forward UID: " << input->GetUniqueID()
                 << " size: " << nodes.size()
                 << " L: " << forwardMin
                 << " M: " << forwardMass
                 << " S: " << forwardScale
                 << " E: " << forwardOffset);

#ifdef USE_BACKWARD
    fitFunction.fForwardFit = false;
    rootMinimizer->SetVariable(0,"log(Mass)", std::log(105.0*unit::MeV), 0.1);
    rootMinimizer->SetVariable(1,"log(Scale)", 0.0, 0.1);
    rootMinimizer->SetVariable(2,"kinOffset", 0.0*unit::MeV, 10*unit::MeV);
    rootMinimizer->Minimize();
    double backwardMin = rootMinimizer->MinValue();
    double backwardMass = std::exp(rootMinimizer->X()[0]);
    double backwardScale = std::exp(rootMinimizer->X()[1]);
    double backwardOffset = std::abs(rootMinimizer->X()[2]);

    CaptNamedLog("TTrackMassFit", 
                 "backward UID: " << input->GetUniqueID()
                 << " size: " << nodes.size()
                 << " L: " << backwardMin
                 << " M: " << backwardMass
                 << " S: " << backwardScale
                 << " E: " << backwardOffset);
#endif

    return input;
}
