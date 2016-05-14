#include <algorithm>

#include "TDistributeCharge.hxx"

#include <CaptGeomId.hxx>
#include <TCaptLog.hxx>

//////////////////////////////////////////////////////////////////////
// TMeasurement
//////////////////////////////////////////////////////////////////////

CP::DistributeCharge::TMeasurement::TMeasurement(const Object& hit,
                                                 double charge)
    : fObject(hit), fCharge(charge) {
}

void CP::DistributeCharge::TMeasurement::Dump(bool dumpLinks) const  {
    CaptLog("TMeasurement(" << std::hex << this << ")"
             << std::dec << " w/ " << fLinks.size() << " links"
            << "  charge: " << GetCharge());
    if (dumpLinks) {
        CP::TCaptLog::IncreaseIndentation();
        for (CP::DistributeCharge::TLinks::const_iterator link=fLinks.begin(); 
             link != fLinks.end(); ++link) {
            (*link)->Dump();
        }
        CP::TCaptLog::DecreaseIndentation();
    }
}

void CP::DistributeCharge::TMeasurement::NormalizeWeights() const {
    double totalWeight = 0;
    double weightOffset = 0.0;
    int throttle = 5;
    do {
        totalWeight = 0;
        for (CP::DistributeCharge::TLinks::const_iterator link = GetLinks().begin();
             link != GetLinks().end(); ++link) {
            double w = (*link)->GetWeight()*(*link)->GetPhysicsWeight();
            if (w<0) w = 0;
            totalWeight += w;
        }
        
        // Make sure that at least some of the weights are positive.  If not,
        // then set some default values.
        if (totalWeight < 1E-6) {
            totalWeight = 0.0;
            for (CP::DistributeCharge::TLinks::const_iterator link
                     = GetLinks().begin();
                 link != GetLinks().end(); ++link) {
                (*link)->SetWeight(1.0/GetLinks().size());
                double w = (*link)->GetWeight()*(*link)->GetPhysicsWeight();
                if (w<0) w = 0;
                totalWeight += w;
            }
        }

        double scaleFactor = (totalWeight-weightOffset)/(1.0-weightOffset);

        weightOffset = 0.0;
        totalWeight = 0.0;
        for (CP::DistributeCharge::TLinks::const_iterator link 
                 = GetLinks().begin();
             link != GetLinks().end(); ++link) {
            double w = (*link)->GetWeight()/scaleFactor;
            if (w<0) w = 0;
            if (w>0.9999) {
                w = 1.0;
                weightOffset += w*(*link)->GetPhysicsWeight();
            }
            totalWeight += w*(*link)->GetPhysicsWeight();
            (*link)->SetWeight(w);
        }
    } while (std::abs(totalWeight-1.0) > 0.001 && 0 <= --throttle);
}

void CP::DistributeCharge::TMeasurement::NormalizePhysicsWeights() const {
    double maxWeight = 0.0;
    for (CP::DistributeCharge::TLinks::const_iterator link = GetLinks().begin();
         link != GetLinks().end(); ++link) {
        double w = (*link)->GetPhysicsWeight();
        if (w<0) w = 0.0;
        maxWeight = std::max(maxWeight,w);
    }
    
    if (maxWeight < 1E-6) {
        for (CP::DistributeCharge::TLinks::const_iterator link = GetLinks().begin();
             link != GetLinks().end(); ++link) {
            (*link)->SetPhysicsWeight(1.0);
        }
        return;
    }

    for (CP::DistributeCharge::TLinks::const_iterator link = GetLinks().begin();
         link != GetLinks().end(); ++link) {
        double w = (*link)->GetPhysicsWeight();
        if (w<0) w = maxWeight;
        (*link)->SetPhysicsWeight(w/maxWeight);
    }

}

double CP::DistributeCharge::TMeasurement::UpdateWeights() const {
    double totalWeight = 0;
    for (CP::DistributeCharge::TLinks::const_iterator link = GetLinks().begin();
         link != GetLinks().end(); ++link) {
        double w = (*link)->GetNewWeight();
        if (w<0) w = 0;
        totalWeight += w;
    }

    if (totalWeight < 1E-6) return 0.0;

    double change = 0;
    double links = 0;
    for (CP::DistributeCharge::TLinks::const_iterator link = GetLinks().begin();
         link != GetLinks().end(); ++link) {
        double w = (*link)->GetNewWeight();
        if (w<0) w = 0;
        w = w/totalWeight;
        change += std::abs(w - (*link)->GetWeight());
        links += 1.0;
        (*link)->SetNewWeight(w);
        (*link)->SetWeight(w);
    }

    if (links>0) change /= links;
    return change;
}

void CP::DistributeCharge::TMeasurement::EliminateLinks(double weightCut) const {
    // There isn't an overlap.
    if (GetLinks().size()<2) return;

    // Get the total weight and the minimum weight.
    double totalWeight = 0;
    double minWeight = 1000.0;
    double maxWeight = 0.0;
    CP::DistributeCharge::TLinks::const_iterator minLink = GetLinks().begin();
    for (CP::DistributeCharge::TLinks::const_iterator link = GetLinks().begin();
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
    for (CP::DistributeCharge::TLinks::const_iterator link 
             = (*minLink)->GetGroup()->GetLinks().begin();
         link !=  (*minLink)->GetGroup()->GetLinks().end();
         ++link) {
        (*link)->SetNewWeight(0.0);
    }
}

void CP::DistributeCharge::TMeasurement::FindLinkWeights() const {
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
    double groupCharge = 0.0;
    for (CP::DistributeCharge::TLinks::const_iterator link = GetLinks().begin();
         link != GetLinks().end(); ++link) {
        double q = (*link)->GetGroup()->GetUniqueCharge(this);
        groupCharge += (*link)->GetGroup()->GetGroupCharge();
        if (q > 0 && (*link)->GetWeight() <= 0) {
            CaptSevere("Zero weight link in group with charge ");
        }
        totalCharge += q;
    }

#ifdef DUMP_DEBUG_INFO
    std::cout << "M " << GetCharge() << " G " << GetCharge() - totalCharge
              << " A " << GetCharge() - groupCharge
              << " " << CP::GeomId::Captain::GetWirePlane(GetObject()->GetGeomId())
              << "-" << CP::GeomId::Captain::GetWireNumber(GetObject()->GetGeomId())
              << std::endl;
#endif

    // Find the new weights for each link in this measurement.  The new weight
    // is the ratio of the unique charge in the measurement group that is
    // linked to to the total charge.
    for (CP::DistributeCharge::TLinks::const_iterator link = GetLinks().begin();
         link != GetLinks().end(); ++link) {
#ifdef DUMP_DEBUG_INFO
        std::cout << "     ";
#endif
        if ((*link)->GetWeight() < 1E-6) {
            (*link)->SetNewWeight(0.0);
        }
        if (totalCharge>1E-9) {
            double w = (*link)->GetWeight();
            double q = (*link)->GetCharge();
            double gq = (*link)->GetGroup()->GetUniqueCharge(this);
            for (CP::DistributeCharge::TLinks::const_iterator lnk
                     = (*link)->GetGroup()->GetLinks().begin();
                 lnk != (*link)->GetGroup()->GetLinks().end(); ++lnk) {
#ifdef DUMP_DEBUG_INFO
                std::cout << CP::GeomId::Captain::GetWireNumber(
                    (*lnk)->GetMeasurement()->GetObject()->GetGeomId());
#endif
            }
            if (q > 1E-6) w *= gq/q;
            else w = 0.0;
            w *= GetCharge()/totalCharge;
            if (w > 1.0) w = 1.0;
#ifdef DUMP_DEBUG_INFO
            std::cout << " q " << q
                      << " g " << gq
                      << " w " << w
                      << " o " << (*link)->GetWeight()
                      << std::endl;
#endif
            (*link)->SetNewWeight(w);
        }
        else {
            (*link)->SetNewWeight(1.0/GetLinks().size());
        }
    }
}

///////////////////////////////////////////////////////////////////////////
// TMeasurementGroup
///////////////////////////////////////////////////////////////////////////

CP::DistributeCharge::TMeasurementGroup::TMeasurementGroup(
    CP::TDistributeCharge* owner, TMeasurementGroup::Object& object) 
    : fOwner(owner), fObject(object) {
}

void CP::DistributeCharge::TMeasurementGroup::Dump(bool dumpLinks) const {
    CaptLog("TMeasurementGroup(" << std::hex << this << ")"
             << std::dec << " w/ " << fLinks.size() << " links"
            << "  charge: " << GetGroupCharge());
    if (dumpLinks) {
        CP::TCaptLog::IncreaseIndentation();
        for (CP::DistributeCharge::TLinks::const_iterator link=fLinks.begin(); 
             link != fLinks.end(); ++link) {
            (*link)->Dump();
        }
        CP::TCaptLog::DecreaseIndentation();
    }
}

double CP::DistributeCharge::TMeasurementGroup::GetGroupCharge() const {
    double charge = 0.0;
    double count = 0.0;
    for (CP::DistributeCharge::TLinks::const_iterator link 
             = GetLinks().begin();
         link != GetLinks().end(); ++link) {
        charge += (*link)->GetCharge();
        count += 1.0;
    }
    if (count < 1.0) return 0.0;
    return charge/count;
}

double CP::DistributeCharge::TMeasurementGroup::GetUniqueCharge(
    const CP::DistributeCharge::TMeasurement* cb) const {
    double charge = 0.0;
    double count = 0.0;
    for (CP::DistributeCharge::TLinks::const_iterator link 
             = GetLinks().begin();
         link != GetLinks().end(); ++link) {
        if ((*link)->GetMeasurement() == cb) continue;
        charge += (*link)->GetCharge();
        count += 1.0;
    }
    if (count < 1.0) return 0.0;
    return charge/count;
}

CP::DistributeCharge::TMeasurement* 
CP::DistributeCharge::TMeasurementGroup::AddMeasurement(
    CP::DistributeCharge::TMeasurement::Object& object,
    double charge,
    double physicsWeight) {
    // Find the existing measurement, or create a new one.
    CP::DistributeCharge::TMeasurement* measurement 
        = fOwner->FindMeasurement(object,charge);
    fOwner->CreateLink(this,measurement,physicsWeight);
    return measurement;
}

void CP::DistributeCharge::TLink::Dump() const {
    CaptLog("TLink(" << std::hex << this << ")"
            << std::dec <<std::setprecision(3) << " Q " << GetRawCharge()
            << std::dec <<std::setprecision(3) << " G "
            << GetGroup()->GetGroupCharge()
            << std::dec <<std::setprecision(3) << " U "
            << GetGroup()->GetUniqueCharge(GetMeasurement())
            << std::dec <<std::setprecision(3) << " q " << GetCharge()
            << std::dec <<std::setprecision(3) << " w " << fWeight
            << std::dec <<std::setprecision(3) << " p " << fPhysicsWeight
            << std::hex << " g " << fMeasurementGroup
            << std::hex << " M " << fMeasurement
            << std::dec);
}

/////////////////////////////////////////////////////////////////////
// TDistributeCharge
/////////////////////////////////////////////////////////////////////
    
CP::TDistributeCharge::TDistributeCharge() 
    : fWeightCut(0.1) {}
CP::TDistributeCharge::~TDistributeCharge() {}

CP::DistributeCharge::TMeasurementGroup& 
CP::TDistributeCharge::AddGroup(
    DistributeCharge::TMeasurementGroup::Object& object) {
    fGroups.push_back(CP::DistributeCharge::TMeasurementGroup(this,object));
    return fGroups.back();
}

CP::DistributeCharge::TMeasurement*
CP::TDistributeCharge::FindMeasurement(
    CP::DistributeCharge::TMeasurement::Object& object,
    double charge) {
#ifdef USE_LIST
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
    fMeasurements.push_back(CP::DistributeCharge::TMeasurement(object,charge));
    return &(fMeasurements.back());
#else
    Measurements::iterator m
        = fMeasurements.insert(
            CP::DistributeCharge::TMeasurement(object,charge)).first;
     return const_cast<CP::DistributeCharge::TMeasurement*>(&(*m));
#endif
}

CP::DistributeCharge::TLink* CP::TDistributeCharge::CreateLink(
    CP::DistributeCharge::TMeasurementGroup* group,
    CP::DistributeCharge::TMeasurement* measurement,
    double physicsWeight) {
    fLinks.push_back(CP::DistributeCharge::TLink(group,measurement));
    CP::DistributeCharge::TLink* link = &fLinks.back();
    link->SetPhysicsWeight(physicsWeight);
    measurement->GetLinks().push_back(link);
    group->GetLinks().push_back(link);
    return link;
}

void CP::TDistributeCharge::DumpGroups(bool dumpLinks) const {
    CaptLog("TDistributeCharge(" << std::hex << this << ")  Groups:");
    CP::TCaptLog::IncreaseIndentation();
    for (Groups::const_iterator g = fGroups.begin(); g != fGroups.end(); ++g) {
        g->Dump(dumpLinks);
    }
    CP::TCaptLog::DecreaseIndentation();
}


void CP::TDistributeCharge::Dump(bool dumpLinks) const {
    DumpGroups(dumpLinks);
}

void CP::TDistributeCharge::DumpMeasurements(bool dumpLinks) const {
    CaptLog("TDistributeCharge(" << std::hex << this << ")  Measurements:");
    CP::TCaptLog::IncreaseIndentation();
    for (Measurements::const_iterator m = fMeasurements.begin();
         m != fMeasurements.end(); ++m) {
        m->Dump(dumpLinks);
    }
    CP::TCaptLog::DecreaseIndentation();
}

double CP::TDistributeCharge::Solve(double tolerance, int iterations) {
    CaptInfo("Share charge with tolerance: " << tolerance);

    // Make sure the input weights are normalized.
    for (Measurements::iterator m = fMeasurements.begin();
         m != fMeasurements.end(); ++m) {
        m->NormalizePhysicsWeights();
    }

    for (Measurements::iterator m = fMeasurements.begin();
         m != fMeasurements.end(); ++m) {
        m->NormalizeWeights();
    }

    // Do the relaxation, but limit the total number of iterations.
    double change = 0.0;
    while (0 < iterations--) {
        change = RelaxWeights();
        if (change < tolerance) break;
    }

    return change;
}

double CP::TDistributeCharge::RelaxWeights() {
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
    
    // Make sure the input weights are normalized.
    for (Measurements::iterator m = fMeasurements.begin();
         m != fMeasurements.end(); ++m) {
        m->NormalizeWeights();
    }

    return delta/fMeasurements.size();
}
