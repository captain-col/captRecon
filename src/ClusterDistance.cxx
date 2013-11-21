#include "ClusterDistance.hxx"

#include <TReconCluster.hxx>
#include <THitSelection.hxx>
#include <THit.hxx>
#include <THandle.hxx>
#include <HEPUnits.hxx>

#include <TVector3.h>

// This finds the minimum distance between hits in the two clusters.
double CP::ClusterDistance(const CP::TReconCluster& a, 
                           const CP::TReconCluster& b) {
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

