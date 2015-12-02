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

void CP::ShareCharge::TMeasurement::EliminateLinks(double weightCut) {
    // There isn't an overlap.
    if (GetLinks().size()<2) return;

    // Get the total weight and the minimum weight.
    double totalWeight = 0;
    double minWeight = 1000.0;
    double maxWeight = 0.0;
    CP::ShareCharge::TLinks::iterator minLink = GetLinks().begin();
    for (CP::ShareCharge::TLinks::iterator link = GetLinks().begin();
         link != GetLinks().end(); ++link) {
        double w = (*link)->GetWeight();
        if (w<0) w = 0;
        if (w>0 && w < minWeight) {
            minLink = link;
            minWeight = w;
        }
        maxWeight = std::max(maxWeight, w);
        totalWeight += w;
    }

    double averageWeight = totalWeight/GetLinks().size();

    if (minWeight > weightCut*averageWeight) return;

    // The minimum link needs to be eliminated.  When it's eliminated, the
    // other links to the TMeasurementGroup have their weights set to zero.
    for (CP::ShareCharge::TLinks::const_iterator link 
             = (*minLink)->GetGroup()->GetLinks().begin();
         link !=  (*minLink)->GetGroup()->GetLinks().end();
         ++link) {
        (*link)->SetNewWeight(0.0);
    }
}

void CP::ShareCharge::TMeasurement::FindLinkWeights() {
    // This should never happen.
    if (GetLinks().size()<0) return;

    // There isn't an overlap.
    if (GetLinks().size()==1) {
        GetLinks().front()->SetNewWeight(1.0);
        return;
    }

    // Find the total charge in all the measurement groups, not including the
    // charge in this measurement.
    double totalCharge = 0.0;
    for (CP::ShareCharge::TLinks::iterator link = GetLinks().begin();
         link != GetLinks().end(); ++link) {
        double q = (*link)->GetGroup()->GetUniqueCharge(this);
        if (q > 0 && (*link)->GetWeight() <= 0) {
            CaptSevere("Zero weight link in group with charge ");
        }
        totalCharge += q;
    }
    // Find the new weights for each link in this measurement.  The new weight
    // is the ratio of the unique charge in the measurement group that is
    // linked to to the total charge.
    for (CP::ShareCharge::TLinks::iterator link = GetLinks().begin();
         link != GetLinks().end(); ++link) {
        if ((*link)->GetWeight() < 1E-6) {
            (*link)->SetNewWeight(0.0);
        }
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
    : fOwner(owner), fObject(object) {
}

void CP::ShareCharge::TMeasurementGroup::Dump(bool dumpLinks) const {
    CaptLog("TMeasurementGroup(" << std::hex << this << ")"
             << std::dec << " w/ " << fLinks.size() << " links"
            << "  charge: " << GetGroupCharge());
    if (dumpLinks) {
        CP::TCaptLog::IncreaseIndentation();
        for (CP::ShareCharge::TLinks::const_iterator link = fLinks.begin(); 
             link != fLinks.end(); ++link) {
            (*link)->Dump();
        }
        CP::TCaptLog::DecreaseIndentation();
    }
}

double CP::ShareCharge::TMeasurementGroup::GetGroupCharge() const {
    double charge = 0.0;
    double count = 0.0;
    for (CP::ShareCharge::TLinks::const_iterator link 
             = GetLinks().begin();
         link != GetLinks().end(); ++link) {
        charge += (*link)->GetCharge();
        count += 1.0;
    }
    if (count < 1.0) return 0.0;
    return charge;
}

double CP::ShareCharge::TMeasurementGroup::GetUniqueCharge(
    const CP::ShareCharge::TMeasurement* cb) const {
    double charge = 0.0;
    double count = 0.0;
    for (CP::ShareCharge::TLinks::const_iterator link 
             = GetLinks().begin();
         link != GetLinks().end(); ++link) {
        if ((*link)->GetMeasurement() == cb) continue;
        charge += (*link)->GetCharge();
        count += 1.0;
    }
    if (count < 1.0) return 0.0;
    return charge;
}

CP::ShareCharge::TMeasurement* 
CP::ShareCharge::TMeasurementGroup::AddMeasurement(
    CP::ShareCharge::TMeasurement::Object& object, double charge) {
    // Find the existing measurement, or create a new one.
    CP::ShareCharge::TMeasurement* measurement 
        = fOwner->FindMeasurement(object,charge);
    fOwner->CreateLink(this,measurement);
    return measurement;
}

void CP::ShareCharge::TLink::Dump() const {
    CaptLog("TLink(" << std::hex << this << ")"
            << std::dec <<std::setprecision(3) << " C " << GetCharge()
            << std::dec <<std::setprecision(3) << " U "
            << GetGroup()->GetUniqueCharge(GetMeasurement())
            << std::dec <<std::setprecision(3) << " q " << GetRawCharge()
            << std::dec <<std::setprecision(3) << " w " << fWeight
            << std::dec <<std::setprecision(3) << " p " << fPhysicsWeight
            << std::hex << " g " << fMeasurementGroup
            << std::hex << " M " << fMeasurement
            << std::dec);
}

/////////////////////////////////////////////////////////////////////
// TShareCharge
/////////////////////////////////////////////////////////////////////
    
CP::TShareCharge::TShareCharge() 
    : fWeightCut(0.1) {}
CP::TShareCharge::~TShareCharge() {}

CP::ShareCharge::TMeasurementGroup& 
CP::TShareCharge::AddGroup(ShareCharge::TMeasurementGroup::Object& object) {
    fGroups.push_back(CP::ShareCharge::TMeasurementGroup(this,object));
    return fGroups.back();
}

CP::ShareCharge::TMeasurement*
CP::TShareCharge::FindMeasurement(CP::ShareCharge::TMeasurement::Object& object,
                                  double charge) {
    for (Measurements::iterator m = fMeasurements.begin();
         m != fMeasurements.end(); ++m) {
        if (m->GetObject() == object) {
            if (std::abs(m->GetCharge()-charge) > 1E-5) {
                CaptError("Charge mismatch");
                throw;
            }
            return &(*m);
        }
    }
    fMeasurements.push_back(CP::ShareCharge::TMeasurement(object,charge));
    return &(fMeasurements.back());
}

CP::ShareCharge::TLink* CP::TShareCharge::CreateLink(
    CP::ShareCharge::TMeasurementGroup* group,
    CP::ShareCharge::TMeasurement* measurement) {
    fLinks.push_back(CP::ShareCharge::TLink(group,measurement));
    CP::ShareCharge::TLink* link = &fLinks.back();
    measurement->GetLinks().push_back(link);
    group->GetLinks().push_back(link);
    return link;
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
    CaptLog("Share charge with tolerance: " << tolerance);

    // Do the relaxation, but limit the total number of iterations.
    double change = 0.0;
    while (0 < iterations--) {
        change = RelaxWeights();
        if (change < tolerance) break;
    }

    CaptLog("Share charge final change: " << change);

    return change;
}

double CP::TShareCharge::RelaxWeights() {
    // Make sure the input weights are normalized.
    for (Measurements::iterator m = fMeasurements.begin();
         m != fMeasurements.end(); ++m) {
        m->NormalizeWeights();
    }

    // Check if any links need to be eliminated.
    for (Measurements::iterator m = fMeasurements.begin();
         m != fMeasurements.end(); ++m) {
        m->EliminateLinks(fWeightCut);
    }
    
    // Update the weights with the changes.
    for (Measurements::iterator m = fMeasurements.begin();
         m != fMeasurements.end(); ++m) {
        m->UpdateWeights();
    }
    
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
    
    return delta/fMeasurements.size();
}
