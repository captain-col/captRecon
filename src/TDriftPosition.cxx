#include <TDriftPosition.hxx>

#include <TVector3.h>

CP::TDriftPosition::TDriftPosition(double z, double v)
    : fZPlane(z), fDriftVelocity(v) {}

TLorentzVector 
CP::TDriftPosition::operator() (const CP::THandle<CP::THit>& hit) {
    TVector3 p(hit->GetPosition());
    double t(hit->GetTime());
    double deltaZ = fZPlane-p.Z();
    double deltaT = deltaZ/fDriftVelocity;
    return TLorentzVector(p.X(), p.Y(), fZPlane, t + deltaT);
}

TLorentzVector 
CP::TDriftPosition::operator() (const CP::THandle<CP::THit>& hit, double tNew) {
    TVector3 p(hit->GetPosition());
    double t(hit->GetTime());
    double deltaT = tNew - t;
    double deltaZ = deltaT * fDriftVelocity;
    return TLorentzVector(p.X(), p.Y(), p.Z() + deltaZ, tNew);
}

double CP::TDriftPosition::GetTime(const CP::THandle<CP::THit>& hit) {
    TVector3 p(hit->GetPosition());
    double t(hit->GetTime());
    double deltaZ = fZPlane-p.Z();
    double deltaT = deltaZ/fDriftVelocity;
    return t + deltaT;
}
