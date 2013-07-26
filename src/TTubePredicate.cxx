#include "TTubePredicate.hxx"

#include <cmath>

CP::TTubePredicate::TTubePredicate(const TVector3& end1, const TVector3& end2,
                                   double vertRad, double horiRad, 
                                   double horiScale) {
    fEnd = end1;
    fDir = end2-fEnd;
    fLength = fDir.Mag();
    fDir = fDir.Unit();
    
    // Define the horizontal direction for the tube.
    fHoriDir = fDir.Cross(TVector3(0,0,1));
    double sinZ = fHoriDir.Mag();
    if (sinZ < 0.1) fHoriDir = TVector3(1,0,0);
    else fHoriDir = fHoriDir.Unit();
    
    // Define the vertical direction for the tube.  If the tube is vertical,
    // then this is just the Y axis.
    fVertDir = fDir.Cross(fHoriDir);

    double cosZ = 1.0-sinZ*sinZ;
    if (cosZ < 1.0) cosZ = std::sqrt(cosZ);
    else cosZ = 1.0;

    // Find the horizontal radius taking into account the width of confusion
    // for a horizontal tube.
    double extraWidth = fLength;
    if (cosZ > 0.001) extraWidth = std::min(extraWidth,horiScale/cosZ);
    fHoriRad = std::max(horiRad, extraWidth);
    
    // Find the vertical radius. 
    fVertRad = std::min(vertRad/std::max(sinZ,0.001),
                        horiRad/std::max(cosZ,0.001));

}

bool CP::TTubePredicate::operator () (const TVector3& pnt) const {
    TVector3 delt = pnt - fEnd;

    // Check that the point is between the end points.
    double len = delt * fDir;
    if (len < 0.0) return false;
    if (len > fLength) return false;

    len = std::abs(delt*fVertDir);
    if (len > fVertRad) return false;

    len = std::abs(delt*fHoriDir);
    if (len > fHoriRad) return false;

    return true;
}
