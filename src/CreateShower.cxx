#include "CreateShower.hxx"

#include <TCaptLog.hxx>
#include <TReconCluster.hxx>
#include <TReconShower.hxx>
#include <THandle.hxx>
#include <ostreamTVector3.hxx>

/// Create an empty unfit shower using a cluster.
CP::THandle<CP::TReconShower> 
CP::CreateShower(const char* name, CP::THandle<CP::TReconCluster> cluster) {
    TVector3 approxDir = cluster->GetLongAxis();
    return CP::CreateShower(name,
                            cluster->GetHits()->begin(), 
                            cluster->GetHits()->end(),
                            approxDir);
}

/// Determine the shower width based on a cluster.
double CP::CreateShowerWidth(CP::THandle<CP::TReconCluster> cluster, 
                             const TVector3& approxDir) {
    TVector3 dir = approxDir;
    dir = dir.Unit();

    TVector3 longAxis = cluster->GetLongAxis();
    double longW = (longAxis - (dir*longAxis)*dir).Mag();
    TVector3 majorAxis = cluster->GetMajorAxis();
    double majorW = (majorAxis - (dir*majorAxis)*dir).Mag();
    TVector3 minorAxis = cluster->GetMinorAxis();
    double minorW = (minorAxis - (dir*minorAxis)*dir).Mag();

    // The shower width is the middle of longW, majorW, and
    // minorW
    double width = 0.0;
    if (longW <= majorW && minorW <= longW) width = longW;
    else if (longW <= minorW && majorW <= longW) width = longW;
    else if (majorW <= longW && minorW <= majorW) width = majorW;
    else if (majorW <= minorW && longW <= majorW) width = majorW;
    else if (minorW <= longW && majorW <= minorW) width = minorW;
    else if (minorW <= majorW && longW <= minorW) width = minorW;
    else {
        CaptError("No possible width");
        return 3.0*unit::mm;
    }
    return width;
}

