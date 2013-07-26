#ifndef TTubePredicate_hxx_seen
#define TTubePredicate_hxx_seen

#include <TVector3.h>
#include <HEPUnits.hxx>

namespace CP {
    class TTubePredicate;
};

/// A class to determine if a point is inside a tube that is defined by the
/// end points.  The definition of the tube is optimized for a LAr like
/// detector where the Z axis is along the drift direction and has the best
/// accuracy.  The horizontal axis is assumed to have a fairly coarse
/// measurement accuracy, and (importantly) to have confusion between points
/// at the same Z.  It doesn't strongly affect the internal math, but
/// conceptually, the hits are assumed to be composed of three crossed UVX
/// wires. The end points are at the center of the tube end caps.  This takes
/// several parameters.  The end points are not ordered.  The vertRad gives
/// the "radius" of the tube in the Z direction.  The horiRad gives a fixed
/// radius in the XY plane.  It determines the minimum horizontal radius
/// independent of the tube direction.  Since there is "hit confusion" that
/// increases as the tube becomes more horizontal, the horiScale is used to
/// control the amount of confusion.  It defines the "Z" distance over which
/// XY hits will be confused (typically about 1 mm), and is used to define a
/// direction dependent horizontal radius within which points are included.
/// There are default values for vertRad, horiRad and horiScale that are
/// approximately right assuming the wires have a 3 mm spacing and the drift
/// distance is less than 2 meters.
class CP::TTubePredicate {
public:
    /// See the class documentation for details.
    TTubePredicate(const TVector3& end1, const TVector3& end2,
                   double vertRad = 3*unit::mm, 
                   double horiRad = 5*unit::mm, 
                   double horiScale = 2*unit::mm);

    /// The predicate to determine if a point is inside a tube.
    bool operator () (const TVector3& pnt) const;
    
private:

    /// One end of the tube.
    TVector3 fEnd;

    /// The direction of the tube.
    TVector3 fDir;
    
    /// The length of the tube.
    double fLength;

    /// The "vertical" normal to the tube.  It is coplanar with the Z axis.
    TVector3 fVertDir;

    /// The radius along the vertical direction.
    double fVertRad;

    /// The horizontal normal to the tube.  
    TVector3 fHoriDir;

    ///  The radius along the horizontal directon.  This is not normalized
    /// since it also contains the radius information.  It combines
    /// information from the end points, the horiRad, and the horiScale.
    double fHoriRad;
};
#endif
