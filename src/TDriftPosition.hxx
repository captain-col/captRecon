#ifndef TDriftPosition_hxx_seen
#define TDriftPosition_hxx_seen

#include <THit.hxx>
#include <HEPUnits.hxx>

#include <TLorentzVector.h>

namespace CP {
    class TDriftPosition;
};

/// A utility class to handle the effects of the time drift in the LAr.  
class CP::TDriftPosition {
public:
    TDriftPosition();
    
    /// Return the x,y,z,t for a hit after drifting to a time, t.  The time is
    /// "absolute" (relative to the trigger), not a delta time.  The total
    /// drift time will be (t-hit->GetTime()).
    TLorentzVector GetPosition(const CP::THit& hit, double t);

    /// Return the corrected time for a hit after drifting to a Z-plane.
    double GetTime(const CP::THit& hit, double z = 0.0);

private:
        
    /// The drift velocity.  This is set by captRecon.driftVelocity in the
    /// parameters file.
    double fDriftVelocity;
};
#endif
