///
/// An unbelievably ugly hack to test private methods.  This needs to be done
/// before the class being tested,
///
#define private public
#define protected public
#include <TEnergyLoss.hxx>
#undef private
#undef protected

#include <HEPUnits.hxx>
#include <TCaptLog.hxx>

#include <tut.h>

namespace tut {
    struct baseEnergyLoss {
        baseEnergyLoss() {
            // Run before each test.
        }
        ~baseEnergyLoss() {
            // Run after each test.
        }
    };

    // Declare the test
    typedef test_group<baseEnergyLoss>::object testEnergyLoss;
    test_group<baseEnergyLoss> groupEnergyLoss("TEnergyLoss");

    // Test the declaration.
    template<> template<> void testEnergyLoss::test<1> () {
        CP::TEnergyLoss eloss;

        ensure_equals("Default material name", 
                      eloss.GetMaterial(), 
                      "captain");
        ensure_tolerance("captain density", 
                         eloss.GetDensity(), 
                         1.396*unit::gram/unit::cm3, 0.0001);
        ensure_tolerance("captain Z/A ratio",
                         eloss.GetZA()*unit::gram/unit::mole,
                         18.0/39.95, 0.0001);
        ensure_tolerance("captain excitation energy", 
                         eloss.GetExcitationEnergy(),
                         10.2*unit::eV, 0.0001);
        ensure_tolerance("captain plasma energy", 
                         eloss.GetPlasmaEnergy(),
                         22.84*unit::eV, 0.0001);
    }

    // Check the most probable energy loss.
    template<> template<> void testEnergyLoss::test<2> () {
        CP::TEnergyLoss eloss;
        ensure_tolerance("Most probable energy loss",
                         eloss.GetMostProbable(std::log(3.0),1.0*unit::cm),
                         2.117*unit::MeV, 0.001);
    }

};

// Local Variables:
// mode:c++
// c-basic-offset:4
// End:
