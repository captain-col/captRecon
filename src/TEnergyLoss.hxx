#ifndef TEnergyLoss_hxx_seen
#define TEnergyLoss_hxx_seen

#include <string>
#include <cmath>

namespace CP {
    class TEnergyLoss;
};

/// Encapsulate the standard approximate energy loss formulas for particles in
/// materials.  This implements the approximate energy loss formulas for
/// particles in materials for use in a reconstruction, but don't depend on it
/// for a precise calculation.  It provides the usual Bethe-Bloch formula, as
/// well as the most probable energy loss and the Landau scale factor.  The
/// effect of detector resolution is included in the Landau scale factor by
/// adding the resolution in quadrature with the scale (at best a crude
/// approximation).  The material properties can be set when the object is
/// constructed, and it knows about a few standard materials (to make it
/// easier to use).
///
/// This uses HEP Units, with the convention that zero mass particles are
/// electrons.  When the particle mass is above fHeavyThreshold, the heavy
/// particle formulas are used.  When the mass is below fHeavyThreshold, the
/// energy loss will be interpolated between the heavy particle energy loss
/// and the electron energy loss.  The energy loss is calculated using the
/// absolute value of the mass.
class CP::TEnergyLoss {
public:
    
    /// Construct a new named energy loss table.  The name is only for
    /// documentation (it's value doesn't matter).  If the material string
    /// matchs one of the know materials then the class is initialized for
    /// that material, otherwise it's just documentation and the user needs to
    /// set the material properties.  If the material value is NULL, then the
    /// default material is used (LAr in CAPTAIN).  If the cutoff energy is
    /// provided, then it is used when calculating the average restricted
    /// energy loss.  Otherwise, the maximum delta-ray energy is used as the
    /// cutoff energy.
    TEnergyLoss(const char* material=NULL, double cutoff=-1.0);
    virtual ~TEnergyLoss();

    /// Get the material name.
    const std::string& GetMaterial() const {return fMaterialName;}

    /// Calculate the average restricted energy loss given the kinetic energy
    /// and mass.  The returned energy loss has units of [energy]/[length]
    /// (i.e. MeV/mm).
    double GetAverage(double kinEnergy, double mass);

    /// Calculate the Bethe-Blohc restricted energy loss.  The returned energy
    /// loss has units of [energy]/[length] (i.e. MeV/mm).  This uses the
    /// "low-energy" approximation for the maximum energy transfer By low
    /// energy that means that 2*gamma*m_e < M where M is the mass of the
    /// particle.  For a muon, that's OK up to about 10 GeV.  It's a 6% error
    /// at 100 GeV (not very large compared to the other approximations in this
    /// class).
    double GetBetheBloch(double logGamma);

    /// Calculate the average ionization loss for an electron.  The returned
    /// energy loss has units of [energy]/[length] (i.e. MeV/mm).  This makes
    /// about 10% difference for the minimum ionizing dEdX (i.e. around gamma
    /// of 3).
    double GetElectronLoss(double logGamma);

    /// Calculate the most probable energy loss when passing through thick
    /// layer of material.  The thickness is given in [length] (i.e. mm).
    /// Keep in mind the caveats that this starts to break down for very thin
    /// radiators (thickness < ~0.1 gram/cm^2 or about 1 mm in water).  The
    /// most probable value for very thin radiators is slightly overestimated.
    double GetMostProbable(double kinEnergy, double mass, 
                           double thickness) const;
    
    /// Calculate the most probable energy loss when passing through a thick
    /// layer of material. 
    double GetMostProbable(double logGamma, double thickness) const;

    /// Calculate the scale factor for the energy fluctuations passing through
    /// a material.  This is basically the FWHM/4 of the Landau
    /// distribution. The thickness is given in [length] (i.e. mm). Keep in
    /// mind the caveats that this starts to break down for very thin
    /// radiators (thickness < ~0.1 gram/cm^2 or about 1 mm in water).  The
    /// fluctuations for very thin radiators is slightly underestimated.
    double GetScaleFactor(double kinEnergy, double mass, double thickness) const;

    /// Calculate the scale factor for the energy fluctuations passing through
    /// a material.  This is basically the FWHM/4 of the Landau distribution.
    /// The thickness is given in [length] (i.e. mm).  Keep in mind the
    /// caveats that this starts to break down for very thin radiators
    /// (thickness < ~0.1 gram/cm^2 or about 1 mm in water).  The fluctuations
    /// for very thin radiators is slightly underestimated.
    double GetScaleFactor(double logGamma, double thickness) const;

    /// Get the PDF value for observing "deposit" energy for a particle.  This
    /// returns the PDF value of observing the deposited energy given a
    /// particle with "kinEnergy" and "mass" passing through a "thickness" of
    /// material.  This takes into account the detector energy resolution.
    double GetDepositPDF(double deposit, double kinEnergy, double mass, 
                         double thickness) const;

    /// Get the PDF value for observing "deposit" energy for a particle.  This
    /// returns the PDF value for observing the deposited energy given a
    /// particle with "logGamma" passing through a "thickness" of material.
    /// This takes into account the detector energy resolution (but the
    /// resolution convolution is approximate).
    double GetDepositPDF(double deposit, double logGamma, 
                         double thickness) const;

    /// Set the cut off energy used in the restricted energy calculation
    void SetCutoffEnergy(double cutoff) {fCutoffEnergy = cutoff;}

    /// Get the current cutoff energy used in the restricted energy loss
    /// calculation.  
    double GetCutoffEnergy() const {return fCutoffEnergy;}

    /// Set material properties.  The fraction is the number fraction of the
    /// atoms in the material.  The Z and A are the usual definitions.  The
    /// excitationEnergy and densityCorrection are as defined in the
    /// BetheBloch formula.  If the excitationEnergy and densityCorrection are
    /// zero, the "PDG" estimates are used.  The default excitation energy is
    /// the Barkas estimate (16Z^(0.9) eV for Z>1 and 19.2 eV for Z=1 [i.e. H2
    /// gas]).
    void SetMaterial(double frac, double z, double a, 
                     double excitationEnergy=0);

    /// Set the density of the material in units of [mass]/[volume].
    void SetDensity(double density) {fDensity = density;}

    /// Get the density of the material.
    double GetDensity() const {return fDensity;}

    /// Get the current Z/A ratio.
    double GetZA() const {return fCurrentZA;}
     
    /// Set the current Z/A ratio.  This overrides the calculated value, and
    /// sets the material fraction to 1.0.
    void SetZA(double chargeMassRatio);

    /// Get the current excitation energy.
    double GetExcitationEnergy() const {return fCurrentExcitationEnergy;}

    /// Set the current excitation energy.  This overrides the estimated value
    /// and sets the material fraction to 1.0.
    void SetExcitationEnergy(double excitationEnergy);

    /// Set the energy resolution parameters.  The percentage resolution is
    /// a+b/sqrt(E)+c/E where E is the deposited energy.
    void SetEnergyResolution(double a, double b, double c) {
        fConstantTerm = a;
        fSqrtTerm = b;
        fLinearTerm = c;
    }

    /// Get the energy resolution.  The energy resolution is returned in units
    /// of energy.  Internally, it's added in quadrature with the Landau scale
    /// factor to determine the scale factor of the "resolution adjusted
    /// Landau".
    double GetEnergyResolution(double energyDeposit) const {
        double res = fConstantTerm 
            + fSqrtTerm/std::sqrt(energyDeposit)
            + fLinearTerm/energyDeposit;
        return res*energyDeposit;
    }
            
    /// Get the current density correction assuming very high energy.
    double GetDensityCorrection(double gamma) const;

    /// Get the plasma energy for the material.  The plasma energyh is a
    /// function of the density, and average Z/A.
    double GetPlasmaEnergy() const;

private:
    /// The current material name
    std::string fMaterialName;

    /// The cutoff energy
    double fCutoffEnergy;

    /// The material density
    double fDensity;

    /// The current fractional weight of the material.  This is only used when
    /// adding new materials.  When adding materials, it's the ratio between
    /// the weights that will matter, but ideally they sum to 1.0.
    double fCurrentWeightSum;

    /// The current average Z/A
    double fCurrentZA;

    /// The current average ionization energy
    double fCurrentExcitationEnergy;

    /// The energy resolution terms
    double fConstantTerm;
    double fSqrtTerm;
    double fLinearTerm;
};
#endif
