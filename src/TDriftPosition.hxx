#ifndef TDriftPosition_hxx_seen
#define TDriftPosition_hxx_seen

#include <THit.hxx>
#include <HEPUnits.hxx>

#include <TLorentzVector.h>

namespace CP {
    class TDriftPosition;
};

class CP::TDriftPosition {
public:
    TDriftPosition(double z=0.0*unit::mm, 
                   double v=1.6*unit::mm/unit::us);
    
    /// Return the x,y,z,t for a hit that has been drifted to a particular Z
    /// plane.  The Z plane will be the plane set in the constructor.  The
    /// time of the hit will be adjusted based on the drift velocity.
    TLorentzVector operator() (const CP::THandle<CP::THit>& hit);

    /// Drift a hit position by the requested amount of time.
    TLorentzVector operator() (const CP::THandle<CP::THit>& hit, double t);

    /// Return the corrected time for a hit after drifting to the plane.
    double GetTime(const CP::THandle<CP::THit>& hit);

private:
        
    /// The position that the hit is drifted to.
    double fZPlane;

    // The drift velocity
    double fDriftVelocity;
};
#endif
