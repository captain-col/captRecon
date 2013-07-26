#include <TClusterMomentsFit.hxx>

#include <HEPUnits.hxx>

#include <TVectorF.h>
#include <TMatrixF.h>

CP::TClusterMomentsFit::TClusterMomentsFit(const CP::TReconCluster& cluster) {
    fPosition = cluster.GetPosition();

    // Find the orientation of the cluster from the eigen vectors.  The first
    // eigen vector (i.e. the long axis) is the "direction" of the cluster.
    TVectorF eigenValues;
    TMatrixF eigenVectors(cluster.GetMoments().EigenVectors(eigenValues));
    fDirection.SetXYZ(eigenVectors(0,0), eigenVectors(1,0), eigenVectors(2,0));
    fDirection = fDirection.Unit();

    const double epsilon = 1E-6;
    if (fDirection.X() < -epsilon) {
        fDirection = -fDirection;
    }
    else if (fDirection.X() < epsilon && fDirection.Y() < -epsilon) {
        fDirection = -fDirection;
    }
    else if (fDirection.Y() < epsilon && fDirection.Z() < 0.0) {
        fDirection = -fDirection;
    }

    // Find the end points.
    double maxZ = 0.0;
    double minZ = 0.0;
    for (CP::THitSelection::const_iterator h = cluster.GetHits()->begin();
         h != cluster.GetHits()->end(); ++h) {
        TVector3 diff = (*h)->GetPosition() - fPosition.Vect();
        double dist = diff*fDirection;
        if (maxZ < dist) maxZ = dist;
        if (minZ > dist) minZ = dist;
    }
    fEnds.first = fPosition 
        + minZ * TLorentzVector(fDirection,unit::ns/(30*unit::cm));
    fEnds.second = fPosition 
        + maxZ * TLorentzVector(fDirection,unit::ns/(30*unit::cm));
}
