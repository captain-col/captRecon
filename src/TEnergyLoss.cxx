#include "TEnergyLoss.hxx"

#include <TCaptLog.hxx>
#include <HEPUnits.hxx>

#include <TMath.h>

#include <cmath>

CP::TEnergyLoss::TEnergyLoss(const char* material, double cutoff) {

    if (material) {
        fMaterialName = std::string(material);
    }
    else {
        fMaterialName = "captain";
    }

    // Setup for the material averages.
    fCurrentWeightSum = 0.0;
    fCurrentZA = 1.0/(2.0*unit::gram/unit::mole);
    fCurrentExcitationEnergy = 10.0*unit::eV;

    // By default, the detector resolution is not perfect.
    fConstantTerm = 0.0;
    fSqrtTerm = 0.0;
    fLinearTerm = 0.0;

    // Set the cutoff energy for the detector.  This is the energy above
    // which delta-rays are resolved as independent particles.
    fCutoffEnergy = cutoff;

    // Set the default density
    fDensity = 1.0*unit::gram/unit::cm3;

    if (fMaterialName == "LAr") {
        SetMaterial(1.0, 18, 39.95*unit::gram/unit::mole, 10.2*unit::eV);
        SetDensity(1.396*unit::gram/unit::cm3);
    }
    else if (fMaterialName == "captain") {
        SetMaterial(1.0, 18, 39.95*unit::gram/unit::mole, 10.2*unit::eV);
        SetDensity(1.396*unit::gram/unit::cm3);
        SetEnergyResolution(0.05,0.0,0.0);
    }
    else if (fMaterialName == "scint") {
        SetMaterial(2, 1, 1.0*unit::gram/unit::mole, 19.2*unit::eV);
        SetMaterial(1, 6, 12*unit::gram/unit::mole, 12*unit::eV);
        SetDensity(1.4*unit::gram/unit::cm3);
    }
}

CP::TEnergyLoss::~TEnergyLoss() {}

void CP::TEnergyLoss::SetMaterial(double frac, double z, double a, 
                                  double excitationEnergy) {
    if (excitationEnergy < 0.1*unit::eV) {
        if (z>1.0) excitationEnergy = 16*std::pow(z,0.9)*unit::eV;
        else excitationEnergy = 19.2*unit::eV;
    }

    if (a < 1000) {
        CaptError("Unreasonable atomic mass: It should be in HEPUnits");
        CaptError("Value is A = " << a);
        CaptError("Corresponds to A = " << a*unit::mole/unit::gram <<" g/mole");
        CaptError("Probably should use  \"A*unit::gram/unit::mole\"");
    }

    double za = fCurrentWeightSum*fCurrentZA + frac*(z/a);
    double ee = fCurrentWeightSum*fCurrentExcitationEnergy
        + frac*excitationEnergy;
    fCurrentWeightSum += frac;
    fCurrentZA = za/fCurrentWeightSum;
    fCurrentExcitationEnergy = ee/fCurrentWeightSum;

}

void CP::TEnergyLoss::SetZA(double za) {
    fCurrentWeightSum = 1.0;
    fCurrentZA = za;
};

void CP::TEnergyLoss::SetExcitationEnergy(double ee) {
    fCurrentWeightSum = 1.0;
    fCurrentExcitationEnergy = ee;
};

double CP::TEnergyLoss::GetPlasmaEnergy() const {
    return 28.8*unit::eV*std::sqrt(GetZA()*(unit::g/unit::mole)
                                   *GetDensity()/(unit::g/unit::cm3));
}

double CP::TEnergyLoss::GetDensityCorrection(double gamma) const {
    if (gamma < 1.001) gamma = 1.001;
    double betaGamma = gamma*std::sqrt(1.0-1.0/(gamma*gamma));
    double hPlasma = GetPlasmaEnergy();
    double halfDelta = std::log(hPlasma/GetExcitationEnergy())
        + std::log(betaGamma)
        - 0.5;
    if (halfDelta < 0.0) halfDelta = 0.0;
    return 2.0*halfDelta;
}

double CP::TEnergyLoss::GetMostProbable(double kinEnergy, double mass,
                                        double thickness) const {
    mass = std::max(mass, 0.510*unit::MeV);
    mass = std::abs(mass);
    kinEnergy = std::abs(kinEnergy);
    double logGamma = std::log((mass+kinEnergy)/mass);
    return GetMostProbable(logGamma, thickness);
}

double CP::TEnergyLoss::GetMostProbable(double logGamma, double thickness) 
    const {
    double K = 0.307075*unit::MeV*unit::cm2/unit::mole;;
    double me = 0.510999*unit::MeV;
    double jFactor = 0.200;
    double gamma = std::exp(std::abs(logGamma));
    double beta2 = 1.0-1.0/(gamma*gamma);
    double zeta = 0.5*K*GetZA()*thickness*fDensity/beta2;
    double arg = std::log(2.0*me*beta2*gamma*gamma/GetExcitationEnergy());
    arg += std::log(zeta/GetExcitationEnergy());
    arg += jFactor;
    arg -= beta2;
    arg -= GetDensityCorrection(gamma);
    return zeta*arg;
}

double CP::TEnergyLoss::GetScaleFactor(double kinEnergy, double mass,
                                       double thickness) const {
    mass = std::max(mass, 0.510*unit::MeV);
    mass = std::abs(mass);
    kinEnergy = std::abs(kinEnergy);
    double logGamma = std::log((mass+kinEnergy)/mass);
    return GetScaleFactor(logGamma, thickness);
}

double CP::TEnergyLoss::GetScaleFactor(double logGamma, double thickness) 
    const {
    double K = 0.307075*unit::MeV*unit::cm2/unit::mole;;
    double gamma = std::exp(std::abs(logGamma));
    double beta2 = 1.0-1.0/(gamma*gamma);
    double zeta = 0.5*K*GetZA()*thickness*fDensity/beta2;
    return zeta;
}

double CP::TEnergyLoss::GetDepositPDF(double deposit, double logGamma,
                                      double thickness) const {
    double mp = GetMostProbable(logGamma,thickness);
    double scale = GetScaleFactor(logGamma,thickness);
    double land = TMath::Landau(deposit,mp,scale,true);
    double resolution = GetEnergyResolution(deposit);
    double gaus = TMath::Gaus(deposit-mp,resolution);
    return 0.999*land + 0.001*gaus + 1E-40;
}

double CP::TEnergyLoss::GetDepositPDF(double deposit, 
                                      double kinEnergy, double mass,
                                      double thickness) const {
    mass = std::abs(mass);
    mass = std::max(mass, 0.510*unit::MeV);
    kinEnergy = std::abs(kinEnergy);
    double logGamma = std::log((mass+kinEnergy)/mass);
    return GetDepositPDF(deposit,logGamma,thickness);
}
    
double CP::TEnergyLoss::GetBetheBloch(double logGamma) {
    double K = 0.307075*unit::MeV*unit::cm2/unit::mole;;
    double gamma = std::exp(std::abs(logGamma));
    double beta2 = 1.0-1.0/(gamma*gamma);
    double zeta = 0.5*K*GetZA()*fDensity/beta2;
    double me = 0.510999*unit::MeV;
    double wMax = 2.0*me*gamma*gamma*beta2;
    double wCut = std::min(GetCutoffEnergy(), wMax);
    double arg = std::log(2.0*me*beta2*gamma*gamma*wCut);
    arg -= 2.0*std::log(GetExcitationEnergy());
    arg -= beta2*(1.0 + wCut/wMax);
    arg -= GetDensityCorrection(gamma);
    return zeta*arg;
}

double CP::TEnergyLoss::GetElectronLoss(double logGamma) {
    double K = 0.307075*unit::MeV*unit::cm2/unit::mole;;
    double gamma = std::exp(std::abs(logGamma));
    double beta2 = 1.0-1.0/(gamma*gamma);
    double zeta = 0.5*K*GetZA()*fDensity/beta2;
    double me = 0.510999*unit::MeV;
    double wMax = 0.5*me*(gamma-1);
    double wCut = std::min(GetCutoffEnergy(), wMax);
    double arg = std::log(2.0*me*beta2*gamma*gamma*wCut);
    arg -= 2.0*std::log(GetExcitationEnergy());
    arg -= beta2*(1.0 + wCut/wMax);
    arg -= GetDensityCorrection(gamma);
    arg += 0.125*((gamma-1.0)/gamma)*((gamma-1.0)/gamma);
    return zeta*arg;
}

double CP::TEnergyLoss::GetAverage(double kinEnergy, double mass) {
    double me = 0.510999*unit::MeV;
    mass = std::max(std::abs(mass), me);
    kinEnergy = std::abs(kinEnergy);
    double logGamma = std::log((mass+kinEnergy)/mass);

    // Find the Bethe-Bloch dEdX and use it for heavy particles.
    double bb = GetBetheBloch(logGamma);

    // Particles over "cutMass" just use Bethe-Bloch.
    double cutMass = 75*unit::MeV;
    if (mass > cutMass) return bb;

    // We have a "light" particle, so find the dEdX for an electron and
    // interpolate between "normal" bethe-bloch and the electron.
    logGamma = std::log(kinEnergy/me);
    double el = GetElectronLoss(logGamma);

    double w = (mass-me)/(cutMass-me);

    return w*bb + (1.0-w)*el;
}

