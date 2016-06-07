#include "ClusterDistance.hxx"

#include <TReconCluster.hxx>
#include <THitSelection.hxx>
#include <THit.hxx>
#include <THandle.hxx>
#include <HEPUnits.hxx>

#include <TVector3.h>

// This finds the minimum distance between hits in the two clusters.
double CP::ClusterDistance(const CP::TReconBase& a, 
                           const CP::TReconBase& b,
                           const int nDists) {
    CP::THandle<CP::THitSelection> aHits = a.GetHits();
    CP::THandle<CP::THitSelection> bHits = b.GetHits();
    double closeDistances[nDists] = {1000*unit::kilometer,
                                        1000*unit::kilometer,
                                        1000*unit::kilometer,
                                        1000*unit::kilometer,
                                        1000*unit::kilometer};
    for (CP::THitSelection::iterator j = aHits->begin(); 
         j != aHits->end(); ++j) {
        for (CP::THitSelection::iterator k = bHits->begin(); 
             k != bHits->end(); ++k) {
            if (j == k) continue;
            double dist = ((*j)->GetPosition()-(*k)->GetPosition()).Mag();
            for (int i = 0; i<nDists; ++i) {
                if (dist > closeDistances[i]) continue;
                std::swap(dist,closeDistances[i]);
            }
        }
    }
    int count = 0;
    for (int i=0; i<nDists; ++i) {
        if (closeDistances[i]<10*unit::m) count = i;
    }
    return closeDistances[count];
}

// This finds the minimum distance between hits in the two clusters.
double CP::MinimumClusterDistance(const CP::TReconBase& a, 
                                  const CP::TReconBase& b) {
    CP::THandle<CP::THitSelection> aHits = a.GetHits();
    CP::THandle<CP::THitSelection> bHits = b.GetHits();
    double minDist = 1000*unit::kilometer;
    for (CP::THitSelection::iterator j = aHits->begin(); 
         j != aHits->end(); ++j) {
        for (CP::THitSelection::iterator k = bHits->begin(); 
             k != bHits->end(); ++k) {
            if (j == k) continue;
            double dist = ((*j)->GetPosition()-(*k)->GetPosition()).Mag();
            if (dist < minDist) minDist = dist;
        }
    }
    return minDist;
}

double CP::ClusterVicinity(const CP::TReconCluster& a, 
                           const CP::TReconCluster& b) {
    const double overlapSigma = 3.0;
    TVector3 tmp[6];
    tmp[0] = a.GetPosition().Vect() + overlapSigma*a.GetLongAxis();
    tmp[1] = a.GetPosition().Vect() - overlapSigma*a.GetLongAxis();
    tmp[2] = a.GetPosition().Vect() + overlapSigma*a.GetMajorAxis();
    tmp[3] = a.GetPosition().Vect() - overlapSigma*a.GetMajorAxis();
    tmp[4] = a.GetPosition().Vect() + overlapSigma*a.GetMinorAxis();
    tmp[5] = a.GetPosition().Vect() - overlapSigma*a.GetMinorAxis();
    TVector3 aMax(-1E10,-1E10,-1E10);
    TVector3 aMin(1E10,1E10,1E10);
    for (int j=0; j<6; ++j) {
        for (int i=0; i<3; ++i) {
            aMax[i] = std::max(tmp[j][i],aMax[i]);
            aMin[i] = std::min(tmp[j][i],aMin[i]);
        }
    }

    tmp[0] = b.GetPosition().Vect() + overlapSigma*b.GetLongAxis();
    tmp[1] = b.GetPosition().Vect() - overlapSigma*b.GetLongAxis();
    tmp[2] = b.GetPosition().Vect() + overlapSigma*b.GetMajorAxis();
    tmp[3] = b.GetPosition().Vect() - overlapSigma*b.GetMajorAxis();
    tmp[4] = b.GetPosition().Vect() + overlapSigma*b.GetMinorAxis();
    tmp[5] = b.GetPosition().Vect() - overlapSigma*b.GetMinorAxis();
    TVector3 bMax(-1E10,-1E10,-1E10);
    TVector3 bMin(1E10,1E10,1E10);
    for (int j=0; j<6; ++j) {
        for (int i=0; i<3; ++i) {
            bMax[i] = std::max(tmp[j][i],bMax[i]);
            bMin[i] = std::min(tmp[j][i],bMin[i]);
        }
    }

    float sqrDist = 0;
    for(int i = 0;i < 3;i++) {
        if(bMax[i] < aMin[i]) {
            float d = bMax[i] - aMin[i];
            sqrDist += d * d;
        }
        else if(bMin[i] > aMax[i])
        {
            float d = bMin[i] - aMax[i];
            sqrDist += d * d;
        }
    }
    if (sqrDist > 0) return std::sqrt(sqrDist);
    return 0.0;
}
