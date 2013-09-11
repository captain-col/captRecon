#include "TDriftPosition.hxx"

#include <THit.hxx>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TRuntimeParameters.hxx>

CP::TDriftPosition::TDriftPosition() {
    fDriftVelocity
        = CP::TRuntimeParameters::Get().GetParameterD(
            "captRecon.driftVelocity");
}

TLorentzVector 
CP::TDriftPosition::GetPosition(const CP::THit& hit, double tNew) const {
    TVector3 p(hit.GetPosition());
    double deltaT = tNew - hit.GetTime();
    double deltaZ = deltaT * fDriftVelocity;
    return TLorentzVector(p.X(), p.Y(), p.Z() + deltaZ, tNew);
}

double CP::TDriftPosition::GetTime(const CP::THit& hit, double zPlane) const {
    double deltaZ = zPlane - hit.GetPosition().Z();
    double deltaT = deltaZ/fDriftVelocity;
    return hit.GetTime() + deltaT;
}
