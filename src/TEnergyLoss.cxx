#include "TEnergyLoss.hxx"

#include <TCaptLog.hxx>
#include <HEPUnits.hxx>
#include <HEPConstants.hxx>

#include <TMath.h>

#include <cmath>

namespace {
    double LogGamma(double kinEnergy, double mass) {
        mass = std::max(std::abs(mass), unit::electron_mass_c2);
        kinEnergy = std::abs(kinEnergy);
        double totalEnergy = mass + kinEnergy;
        double gamma = totalEnergy/mass;
        double logGamma = std::log(gamma);
        return logGamma;
    }
}

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
        SetEnergyResolution(0.03,0.006,30*unit::keV);
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

    fMinimumIonizationGamma = -1;
}

void CP::TEnergyLoss::SetZA(double za) {
    fCurrentWeightSum = 1.0;
    fCurrentZA = za;
    fMinimumIonizationGamma = -1;
};

void CP::TEnergyLoss::SetExcitationEnergy(double ee) {
    fCurrentWeightSum = 1.0;
    fCurrentExcitationEnergy = ee;
    fMinimumIonizationGamma = -1;
};

double CP::TEnergyLoss::GetPlasmaEnergy() const {
    return 28.8*unit::eV*std::sqrt(GetZA()*(unit::g/unit::mole)
                                   *GetDensity()/(unit::g/unit::cm3));
}

double CP::TEnergyLoss::GetEnergyResolution(double energyDeposit) const {
    double rConst = fConstantTerm*energyDeposit;
    double rSqrt = fSqrtTerm*sqrt(energyDeposit);
    double rLinear = fLinearTerm;
    return std::sqrt(rConst*rConst + rSqrt*rSqrt + rLinear*rLinear);
}

double CP::TEnergyLoss::GetDensityCorrection(double gamma) const {
    if (gamma < 1.001) gamma = 1.001;
    if (gamma > 1000) {
        double betaGamma = gamma*std::sqrt(1.0-1.0/(gamma*gamma));
        double hPlasma = GetPlasmaEnergy();
        double halfDelta = std::log(hPlasma/GetExcitationEnergy())
            + std::log(betaGamma)
            - 0.5;
        if (halfDelta < 0.0) halfDelta = 0.0;
        return 2.0*halfDelta;
    }
    // Delta parameterization from PRB.26.6067 Sternheimer, Seltzer and
    // Berger, "Density effect for the ionization of charged particles in
    // various substances".  This is for Argon.
    double beta = std::sqrt(1.0 - 1/(gamma*gamma));
    double X = std::log10(beta*gamma);
    double a = 0.1902;
    double m = 2.982;
    double X0 = 1.1716;
    double X1 = 4.5;
    double C = 11.948;
    if (X <= X0) {
        return 0.0;
    }
    else if (X < X1) {
        return 4.6052*X + a*std::pow(X1-X,m) - C;
    }
    return 4.6052*X - C;
}

double CP::TEnergyLoss::GetMostProbable(double kinEnergy, double mass,
                                        double thickness,
                                        bool truncated) const {
    double logGamma = LogGamma(kinEnergy,mass);
    return GetMostProbable(logGamma, thickness, truncated);
}

// The most probable energy is implemented from PDG 2011 eqn 27.11
double CP::TEnergyLoss::GetMostProbable(double logGamma, double thickness,
                                        bool truncated) const {
    double K = 0.307075;
    double me = unit::electron_mass_c2/unit::MeV;
    double jFactor = 0.200;
    double gamma = std::exp(std::abs(logGamma));
    if (truncated) {
        FindMIPGamma();
        gamma = std::min(gamma,fMinimumIonizationGamma);
    }
    double beta2 = 1.0-1.0/(gamma*gamma);
    double excitationEnergy = GetExcitationEnergy();
    double density = (thickness*fDensity)/(unit::gram/unit::cm2);
    double za = GetZA()/(unit::mole/unit::gram);
    double zeta = 0.5*K*za*density/beta2;
    double arg = std::log(2.0*me*beta2*gamma*gamma/excitationEnergy);
    arg += std::log(zeta/excitationEnergy);
    arg += jFactor;
    arg += - beta2;
    arg += - GetDensityCorrection(gamma);
    return zeta*arg;
}

double CP::TEnergyLoss::GetMIPGamma() const {
    FindMIPGamma();
    return fMinimumIonizationGamma;
}

void CP::TEnergyLoss::FindMIPGamma() const {
    if (fMinimumIonizationGamma > 0) return;
    double mipGamma = 3.0;
    double delta = 1.0;
    double mip = GetMostProbable(std::log(mipGamma), 1*unit::mm,false);
    while (std::abs(delta) > 0.001) {
        double newMIP = GetMostProbable(std::log(mipGamma+delta),1*unit::mm,
                                        false);
        if (newMIP < mip) {
            mip = newMIP;
            mipGamma = mipGamma+delta;
        }
        else if (delta > 0.0) {
            delta = -delta;
        }
        else {
            delta = -0.5*delta;
        }
    }
    fMinimumIonizationGamma = mipGamma;
}

double CP::TEnergyLoss::GetScaleFactor(double kinEnergy, double mass,
                                       double thickness) const {
    double logGamma = LogGamma(kinEnergy,mass);
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

namespace {
    double LandauGaussian(double deposit, 
                          double mpv, double scale, double sigma) {
        // Copied from:
        //
        //   http://root.cern.ch/root/html/tutorials/fit/langaus.C.html
        //
        // In the Landau distribution (represented by the CERNLIB
        // approximation), the maximum is located at x=-0.22278298 with the
        // location parameter=0.  This shift is corrected within this
        // function, so that the actual maximum is identical to the MP
        // parameter.
        
        // Numeric constants
        const double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
        const double mpshift  = -0.22278298;       // Landau maximum location
        
        // Control constants
        const int np = 100;        // number of convolution steps
        const double sc = 5.0;   // convolution extends to +-sc sigmas
        
        // MP shift correction
        double mpc = mpv - mpshift * scale; 
        
        // Range of convolution integral
        double xlow = deposit - sc * sigma;
        if (xlow < 0.0) xlow = 0.0;
        double xup = deposit + sc * sigma;
        double step = (xup-xlow) / np;
        
        // Convolution integral of Landau and Gaussian by sum
        double sum = 0.0;
        for (int i = 1 ; i <= np/2; ++i) {
            double xx = xlow + (1.0*i - 0.5) * step;
            double land = TMath::Landau(xx,mpc,scale) / scale;
            sum += land * TMath::Gaus(deposit,xx,sigma);
            
            xx = xup - (1.0*i - 0.5) * step;
            land = TMath::Landau(xx,mpc,scale) / scale;
            sum += land * TMath::Gaus(deposit,xx,sigma);
        }
        
        double value = (step * sum * invsq2pi / sigma);
        return value;
    }
}

double CP::TEnergyLoss::GetDepositPDF(double logGamma,
                                      double thickness,
                                      double deposit,
                                      double sigma) const {
    double mpv = GetMostProbable(logGamma,thickness);
    double scale = GetScaleFactor(logGamma,thickness);
    double resolution = GetEnergyResolution(deposit);
    resolution = resolution*resolution + sigma*sigma;
    resolution = std::sqrt(resolution);
    return LandauGaussian(deposit, mpv, scale, resolution) + 1E-40;
}

double CP::TEnergyLoss::GetDepositPDF(double kinEnergy, double mass,
                                      double thickness,
                                      double deposit,
                                      double sigma) const {
    double logGamma = LogGamma(kinEnergy,mass);
    return GetDepositPDF(logGamma,thickness,deposit,sigma);
}
    
double CP::TEnergyLoss::GetBetheBloch(double logGamma) {
    double K = 0.307075*unit::MeV*unit::cm2/unit::mole;;
    double gamma = std::exp(std::abs(logGamma));
    double beta2 = 1.0-1.0/(gamma*gamma);
    double zeta = 0.5*K*GetZA()*fDensity/beta2;
    double me = unit::electron_mass_c2;
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
    double me = unit::electron_mass_c2;
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
    mass = std::abs(mass);
    double logGamma = LogGamma(kinEnergy,mass);

    // Find the Bethe-Bloch dEdX and use it for heavy particles.
    double bb = GetBetheBloch(logGamma);

    // Particles over "cutMass" just use Bethe-Bloch.
    double cutMass = 75*unit::MeV;
    if (mass > cutMass) return bb;

    // We have a "light" particle, so find the dEdX for an electron and
    // interpolate between "normal" bethe-bloch and the electron.
    logGamma = LogGamma(kinEnergy,unit::electron_mass_c2);
    double el = GetElectronLoss(logGamma);

    double w = (mass-unit::electron_mass_c2)/(cutMass-unit::electron_mass_c2);
    if (w < 0.0) w = 0.0;
    if (w > 1.0) w = 1.0;

    return w*bb + (1.0-w)*el;
}

