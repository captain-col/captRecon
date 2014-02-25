#include <algorithm>

#include "TShareCharge.hxx"

#include <TCaptLog.hxx>

//////////////////////////////////////////////////////////////////////
// TMeasurement
//////////////////////////////////////////////////////////////////////

CP::ShareCharge::TMeasurement::TMeasurement(const Object& hit, double charge)
    : fObject(hit), fCharge(charge) {
}

void CP::ShareCharge::TMeasurement::Dump(bool dumpLinks) const  {
    CaptLog("TMeasurement(" << std::hex << this << ")"
             << std::dec << " w/ " << fLinks.size() << " links"
            << "  charge: " << GetCharge());
    if (dumpLinks) {
        CP::TCaptLog::IncreaseIndentation();
        for (CP::ShareCharge::TLinks::const_iterator link = fLinks.begin(); 
             link != fLinks.end(); ++link) {
            (*link)->Dump();
        }
        CP::TCaptLog::DecreaseIndentation();
    }
}

void CP::ShareCharge::TMeasurement::NormalizeWeights() {
    double totalWeight = 0;
    for (CP::ShareCharge::TLinks::iterator link = GetLinks().begin();
         link != GetLinks().end(); ++link) {
        double w = (*link)->GetWeight();
        if (w<0) w = 0;
        totalWeight += w;
    }
    
    if (totalWeight > 1E-6) {
        for (CP::ShareCharge::TLinks::iterator link 
                 = GetLinks().begin();
             link != GetLinks().end(); ++link) {
            double w = (*link)->GetWeight();
            if (w<0) w = 0;
            (*link)->SetWeight(w/totalWeight);
        }
    }
    else {
        totalWeight = GetLinks().size();
        for (CP::ShareCharge::TLinks::iterator link
                 = GetLinks().begin();
             link != GetLinks().end(); ++link) {
            (*link)->SetWeight(1.0/totalWeight);
        }
    }
}

double CP::ShareCharge::TMeasurement::UpdateWeights() {
    double totalWeight = 0;
    for (CP::ShareCharge::TLinks::iterator link = GetLinks().begin();
         link != GetLinks().end(); ++link) {
        double w = (*link)->GetNewWeight();
        if (w<0) w = 0;
        totalWeight += w;
    }

    if (totalWeight < 1E-6) return 0.0;

    double change = 0;
    for (CP::ShareCharge::TLinks::iterator link = GetLinks().begin();
         link != GetLinks().end(); ++link) {
        double w = (*link)->GetNewWeight();
        if (w<0) w = 0;
        w = w/totalWeight;
        change += std::abs(w - (*link)->GetWeight());
        (*link)->SetNewWeight(w);
        (*link)->SetWeight(w);
    }

    return change;
}

void CP::ShareCharge::TMeasurement::FindLinkWeights() {
    // This should never happen.
    if (GetLinks().size()<0) return;

    // There isn't an overlap.
    if (GetLinks().size()==1) {
        GetLinks().front()->SetNewWeight(1.0);
        return;
    }

    // Find the total charge in all the clusters, not including this charge
    // bin.
    double totalCharge = 0.0;
    for (CP::ShareCharge::TLinks::iterator link = GetLinks().begin();
         link != GetLinks().end(); ++link) {
        totalCharge += (*link)->GetGroup()->GetUniqueCharge(this);
    }
    // Find the new weights for each link in this charge bin.
    for (CP::ShareCharge::TLinks::iterator link = GetLinks().begin();
         link != GetLinks().end(); ++link) {
        if (totalCharge>1E-9) {
            double w = (*link)->GetGroup()->GetUniqueCharge(this);
            (*link)->SetNewWeight(w/totalCharge);
        
        }
        else {
            (*link)->SetNewWeight(1.0/GetLinks().size());
        }
    }
}

///////////////////////////////////////////////////////////////////////////
// TMeasurementGroup
///////////////////////////////////////////////////////////////////////////

CP::ShareCharge::TMeasurementGroup::TMeasurementGroup(
    CP::TShareCharge* owner, TMeasurementGroup::Object& object) 
    : fOwner(owner), fObject(object) { }

void CP::ShareCharge::TMeasurementGroup::Dump(bool dumpLinks) const {
    CaptLog("TMeasurementGroup(" << std::hex << this << ")"
             << std::dec << " w/ " << fLinks.size() << " links"
            << "  total charge: " << GetTotalCharge());
    if (dumpLinks) {
        CP::TCaptLog::IncreaseIndentation();
        for (CP::ShareCharge::TLinks::const_iterator link = fLinks.begin(); 
             link != fLinks.end(); ++link) {
            (*link)->Dump();
        }
        CP::TCaptLog::DecreaseIndentation();
    }
}

double CP::ShareCharge::TMeasurementGroup::GetTotalCharge() const {
    double charge = 0;
    for (CP::ShareCharge::TLinks::const_iterator link 
             = GetLinks().begin();
         link != GetLinks().end(); ++link) {
        charge += (*link)->GetCharge();
    }
    return charge;
}

double CP::ShareCharge::TMeasurementGroup::GetUniqueCharge(
    const CP::ShareCharge::TMeasurement* cb) const {
    double charge = 0;
    for (CP::ShareCharge::TLinks::const_iterator link 
             = GetLinks().begin();
         link != GetLinks().end(); ++link) {
        if ((*link)->GetMeasurement() == cb) continue;
        charge += (*link)->GetCharge();
    }
    return charge;
}

CP::ShareCharge::TMeasurement& 
CP::ShareCharge::TMeasurementGroup::AddMeasurement(
    CP::ShareCharge::TMeasurement::Object& object, double charge) {
    // Find the existing measurement, or create a new one.
    CP::ShareCharge::TMeasurement& measurement 
        = fOwner->FindMeasurement(object,charge);
    fOwner->CreateLink(*this,measurement);
    return measurement;
}

void CP::ShareCharge::TLink::Dump() const {
    CaptLog("TLink(" << std::hex << this << ")"
            << std::dec <<std::setprecision(3) << " w " << fWeight
            << std::dec <<std::setprecision(3) << " n " << fNewWeight
            << std::dec <<std::setprecision(3) << " p " << fPhysicsWeight
            << std::dec <<std::setprecision(3) << " q " << GetRawCharge()
            << std::dec <<std::setprecision(3) << " C " << GetCharge()
            << std::hex << " g " << fMeasurementGroup
            << std::hex << " M " << fMeasurement
            << std::dec);
}

/////////////////////////////////////////////////////////////////////
// TShareCharge
/////////////////////////////////////////////////////////////////////
    
CP::ShareCharge::TMeasurementGroup& 
CP::TShareCharge::AddGroup(ShareCharge::TMeasurementGroup::Object& object) {
    fGroups.push_back(CP::ShareCharge::TMeasurementGroup(this,object));
    return fGroups.back();
}

CP::ShareCharge::TMeasurement& 
CP::TShareCharge::FindMeasurement(CP::ShareCharge::TMeasurement::Object& object,
                                  double charge) {
    for (Measurements::iterator m = fMeasurements.begin();
         m != fMeasurements.end(); ++m) {
        if (m->GetObject() == object) {
            if (std::abs(m->GetCharge()-charge) > 1E-5) {
                CaptError("Charge mismatch");
                throw;
            }
            return *m;
        }
    }
    fMeasurements.push_back(CP::ShareCharge::TMeasurement(object,charge));
    return fMeasurements.back();
}

CP::ShareCharge::TLink& CP::TShareCharge::CreateLink(
    CP::ShareCharge::TMeasurementGroup& group,
    CP::ShareCharge::TMeasurement& measurement) {
    fLinks.push_back(CP::ShareCharge::TLink(group,measurement));
    return fLinks.back();
}

void CP::TShareCharge::DumpGroups(bool dumpLinks) const {
    CaptLog("TShareCharge(" << std::hex << this << ")  Groups:");
    CP::TCaptLog::IncreaseIndentation();
    for (Groups::const_iterator g = fGroups.begin(); g != fGroups.end(); ++g) {
        g->Dump(dumpLinks);
    }
    CP::TCaptLog::DecreaseIndentation();
}


void CP::TShareCharge::Dump(bool dumpLinks) const {
    DumpGroups(dumpLinks);
}

void CP::TShareCharge::DumpMeasurements(bool dumpLinks) const {
    CaptLog("TShareCharge(" << std::hex << this << ")  Measurements:");
    CP::TCaptLog::IncreaseIndentation();
    for (Measurements::const_iterator m = fMeasurements.begin();
         m != fMeasurements.end(); ++m) {
        m->Dump(dumpLinks);
    }
    CP::TCaptLog::DecreaseIndentation();
}

double CP::TShareCharge::Solve(double tolerance, int iterations) {

    // Do the relaxation, but limit the total number of iterations.
    double change = 0.0;
    while (0 < iterations--) {
        change = RelaxWeights();
        if (change < tolerance) break;
    }

    return change;
}

double CP::TShareCharge::RelaxWeights() {
    // Make sure the input weights are normalized.
    for (Measurements::iterator m = fMeasurements.begin();
         m != fMeasurements.end(); ++m) {
        m->NormalizeWeights();
    }

    // Do one iteration of relaxation.
    for (Measurements::iterator m = fMeasurements.begin();
         m != fMeasurements.end(); ++m) {
        m->FindLinkWeights();
    }

    // Update the weights with the changes.
    double delta = 0.0;
    for (Measurements::iterator m = fMeasurements.begin();
         m != fMeasurements.end(); ++m) {
        delta += m->UpdateWeights();
    }
    
    return delta;
}
