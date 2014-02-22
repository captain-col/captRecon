#include <algorithm>

#include "TP0DShareCharge.hxx"
#include "TP0DAttenuationCurve.hxx"

#include <THandle.hxx>
#include <TReconBase.hxx>
#include <TReconNode.hxx>
#include <TReconVertex.hxx>
#include <TReconShower.hxx>
#include <TVertexState.hxx>
#include <TReconState.hxx>
#include <THitSelection.hxx>
#include <TND280Log.hxx>
#include <FunctionObject.hxx>

#define DUMP_LINKS

/////////////////////////////////////////////////////////////////////
// TP0DShareCharge
/////////////////////////////////////////////////////////////////////

ND::TP0DShareCharge::TP0DShareCharge()
    : TAlgorithm("TP0DShareCharge") {}

ND::TP0DShareCharge::~TP0DShareCharge() {}

ND::THandle<ND::TAlgorithmResult>
ND::TP0DShareCharge::Process(const ND::TAlgorithmResult& input) {

    // Check to see if there is an existing vertex in the input result, or
    // that we have supplemental hits.
    ND::THandle<ND::TReconObjectContainer> inputObjects
        = input.GetResultsContainer();

    // Check that there were input recon base objects. 
    if (!inputObjects) {
        ND280Severe("Must provide TReconObjectContainer in input");
        return ND::THandle<ND::TAlgorithmResult>(); 
    }

    // Check that at least one vertex exists in the input objects.
    bool ok = false;
    for (ND::TReconObjectContainer::iterator o = inputObjects->begin();
         o != inputObjects->end(); ++o) {
        ND::THandle<ND::TReconVertex> vtx = (*o);
        if (vtx) ok = true;
    }
    if (!ok) {
        ND280Severe("An input vertex must be provided.");
        return ND::THandle<ND::TAlgorithmResult>(); 
    }

    // Make a place to save the algorithm results
    ND::THandle<ND::TAlgorithmResult> result(CreateResult());

    // Get the input hits.  This might be empty
    THitSelection* inputHits(new THitSelection("inputHits"));
    if (input.GetHitSelection()) {
        std::copy(input.GetHitSelection()->begin(), 
                  input.GetHitSelection()->end(), 
                  std::back_inserter(*inputHits));
    }
    // Make sure the hits are unique.
    std::sort(inputHits->begin(), inputHits->end(), ND::P0D::geomIdLessThan());
    THitSelection::iterator newEnd = std::unique(inputHits->begin(), 
                                                 inputHits->end(), 
                                                 ND::P0D::geomIdEquals());
    inputHits->erase(newEnd,inputHits->end());
    result->AddHitSelection(inputHits);

    ND::TReconObjectContainer* finalObjects
        = new ND::TReconObjectContainer("final");
    
    result->AddResultsContainer(finalObjects);

    // Go through the input objects and share energy for all of the vertices.
    for (ND::TReconObjectContainer::const_iterator o = inputObjects->begin();
         o != inputObjects->end(); ++o) {
        ND::THandle<ND::TReconVertex> inputVertex = *o;
        if (!inputVertex) {
            // Just pass this through.
            finalObjects->push_back(*o);
        }

        // Create the output vertex so we can add constituents.
        ND::THandle<ND::TReconVertex> outputVertex(new TReconVertex);
        outputVertex->SetAlgorithmName(GetName());
        outputVertex->AddDetector(ND::TReconBase::kP0D);
        outputVertex->ClearStatus(ND::TReconBase::kSuccess);
        ND::THitSelection *vertexHits(new THitSelection("vertexHits"));
        std::copy(inputVertex->GetHits()->begin(), 
                  inputVertex->GetHits()->end(),
                  std::back_inserter(*vertexHits));
        outputVertex->AddHits(vertexHits);
        
        ND::THandle<ND::TVertexState> vertexState 
            = outputVertex->GetState();
        vertexState->SetPosition(inputVertex->GetPosition().X(),
                                 inputVertex->GetPosition().Y(),
                                 inputVertex->GetPosition().Z(),
                                 inputVertex->GetPosition().T());
        ND::THandle<ND::TVertexState> inputState = inputVertex->GetState();
        for (int i=0; i<vertexState->GetDimensions(); ++i) {
            for (int j=0; j<vertexState->GetDimensions(); ++j) {
                vertexState->SetCovarianceValue(
                    i,j,
                    inputState->GetCovarianceValue(i,j));
            }
        }
        outputVertex->SetStatus(ND::TReconBase::kSuccess);

        // Collect all of the showers associated with this vertex.
        ND::THandle<ND::TReconObjectContainer> 
            inputShowers(new TReconObjectContainer);
        for (ND::TReconObjectContainer::iterator c 
                 = inputVertex->GetConstituents()->begin();
             c != inputVertex->GetConstituents()->end(); ++c) {
            ND::THandle<ND::TReconShower> s = *c;
            if (!s) {
                outputVertex->AddConstituent(*c);
                continue;
            }
            inputShowers->push_back(s);
        }

        // Share the energy between the showers
        ND::THandle<ND::TReconObjectContainer> outputShowers
            = ShareCharge(*inputShowers);

        // Add the showers to the vertex and save it in the output objects.
        for (ND::TReconObjectContainer::iterator s = outputShowers->begin();
             s != outputShowers->end(); ++s) {
            outputVertex->AddConstituent(*s);
        }

        finalObjects->push_back(outputVertex);
    }

    return result;
}

ND::THandle<ND::TReconObjectContainer> 
ND::TP0DShareCharge::ShareCharge(const ND::TReconObjectContainer& input) const {
    ND::THandle<ND::TReconObjectContainer> output(new TReconObjectContainer);
    ND::Share::TClusters clusters;
    ND::Share::TCharges charges;
    ND::Share::TLinks links;

    for (ND::TReconObjectContainer::const_iterator shower = input.begin();
         shower != input.end(); ++shower) {
        // Add a new cluster for the shower.
        clusters.AddCluster(*shower,charges,links);
    }

    // Find the optimal weights.  
    for (int i=0; i<5000; ++i) {
        double delta = RelaxWeights(charges);
        if (delta < 1.0E-5) break;
    }

    // Put the results into the output.
    for (ND::Share::TClusters::iterator cluster = clusters.begin();
         cluster != clusters.end(); ++cluster) {
        ND::THandle<ND::TReconShower> shower = cluster->GetReconObject();
        if (shower) {
            ND::THandle<ND::TReconShower> newShower(new TReconShower(*shower));
            newShower->SetAlgorithmName(GetName());
            ND::TReconNodeContainer& newNodes = newShower->GetNodes();
            for (ND::Share::TCluster::iterator clusterBin = cluster->begin();
                 clusterBin != cluster->end(); ++clusterBin) {
                ND::THandle<ND::TReconNode> oldNode = clusterBin->GetNode();
                ND::THandle<ND::TReconNode> newNode;
                // Find the new node that matchs the old node.
                for (ND::TReconNodeContainer::iterator n = newNodes.begin();
                     n != newNodes.end(); ++n) {
                    if (oldNode->GetObject() == (*n)->GetObject()) {
                        newNode = (*n);
                        break;
                    }
                }
                if (!newNode) {
                    ND280Severe("Missing Node in shower");
                    continue;
                }
                ND::THandle<ND::TShowerState> newState = newNode->GetState();
                ND::THandle<ND::TShowerState> oldState = oldNode->GetState();
                newState->SetEDeposit(clusterBin->GetTotalCharge());
                newState->SetEDepositVariance(oldState->GetEDepositVariance());
            }

            // Find the total energy deposit for the shower.
            double dep = 0;
            double var = 0;
            for (ND::TReconNodeContainer::iterator node = newNodes.begin();
                 node != newNodes.end(); ++node) {
                ND::THandle<ND::TShowerState> nodeState = (*node)->GetState();
                dep += nodeState->GetEDeposit();
                var += nodeState->GetEDepositVariance();
            }
            ND::THandle<ND::TShowerState> newState = newShower->GetState();
            newState->SetEDeposit(dep);
            newState->SetEDepositVariance(var);
            
            output->push_back(newShower);
            continue;
        }
        ND280Severe("Missing Shower");
    }

    // Delete the links between the charges.
    for (ND::Share::TLinks::iterator link = links.begin(); 
         link != links.end(); ++link) {
        delete (*link);
    }

    return output;
}

double ND::TP0DShareCharge::RelaxWeights(ND::Share::TCharges& charges) const {
    // Make sure the input weights are normalized.
    for (ND::Share::TCharges::iterator c = charges.begin();
         c != charges.end(); ++c) {
        c->second.NormalizeWeights();
    }

    // Do one iteration of relaxation.
    for (ND::Share::TCharges::iterator c = charges.begin();
         c != charges.end(); ++c) {
        c->second.FindLinkWeights();
    }

    // Update the weights with the changes.
    double delta = 0.0;
    for (ND::Share::TCharges::iterator c = charges.begin();
         c != charges.end(); ++c) {
        delta += c->second.UpdateWeights();
    }
    
    return delta;
}

//////////////////////////////////////////////////////////////////////
// TChargeBin
//////////////////////////////////////////////////////////////////////

ND::Share::TChargeBin::TChargeBin(const ND::THandle<ND::THit>& hit)
    : fHit(hit) {
    fClusterLinks.clear();
}

void ND::Share::TChargeBin::Dump() const  {
    ND280Log("TChargeBin(" << std::hex << this << ")"
             << std::dec << " w/ " << fClusterLinks.size() << " links");
#ifdef DUMP_LINKS
    ND::TND280Log::IncreaseIndentation();
    for (ND::Share::TLinks::const_iterator link = fClusterLinks.begin(); 
         link != fClusterLinks.end(); ++link) {
        (*link)->Dump();
    }
    ND::TND280Log::DecreaseIndentation();
#endif
}

void ND::Share::TChargeBin::NormalizeWeights() {

    double totalWeight = 0;
    for (ND::Share::TLinks::iterator clusterLink = GetClusterLinks().begin();
         clusterLink != GetClusterLinks().end(); ++clusterLink) {
        double w = (*clusterLink)->GetWeight();
        if (w<0) w = 0;
        totalWeight += w;
    }
    
    if (totalWeight > 1E-6) {
        for (ND::Share::TLinks::iterator clusterLink 
                 = GetClusterLinks().begin();
             clusterLink != GetClusterLinks().end(); ++clusterLink) {
            double w = (*clusterLink)->GetWeight();
            if (w<0) w = 0;
            (*clusterLink)->SetWeight(w/totalWeight);
        }
    }
    else {
        totalWeight = GetClusterLinks().size();
        for (ND::Share::TLinks::iterator clusterLink
                 = GetClusterLinks().begin();
             clusterLink != GetClusterLinks().end(); ++clusterLink) {
            (*clusterLink)->SetWeight(1.0/totalWeight);
        }
    }
}

double ND::Share::TChargeBin::UpdateWeights() {
    double totalWeight = 0;
    for (ND::Share::TLinks::iterator clusterLink = GetClusterLinks().begin();
         clusterLink != GetClusterLinks().end(); ++clusterLink) {
        double w = (*clusterLink)->GetNewWeight();
        if (w<0) w = 0;
        totalWeight += w;
    }

    if (totalWeight < 1E-6) return 0.0;

    double change = 0;
    for (ND::Share::TLinks::iterator clusterLink = GetClusterLinks().begin();
         clusterLink != GetClusterLinks().end(); ++clusterLink) {
        double w = (*clusterLink)->GetNewWeight();
        if (w<0) w = 0;
        w = w/totalWeight;
        change += std::abs(w - (*clusterLink)->GetWeight());
        (*clusterLink)->SetNewWeight(w);
        (*clusterLink)->SetWeight(w);
    }

    return change;
}

void ND::Share::TChargeBin::FindLinkWeights() {
    // This should never happen.
    if (GetClusterLinks().size()<0) return;
    // There isn't an overlap.
    if (GetClusterLinks().size()==1) {
        GetClusterLinks().front()->SetNewWeight(1.0);
        return;
    }
    // Find the total charge in all the clusters, not including this charge bin.
    double totalCharge = 0.0;
    for (ND::Share::TLinks::iterator clusterLink = GetClusterLinks().begin();
         clusterLink != GetClusterLinks().end(); ++clusterLink) {
        totalCharge += (*clusterLink)->GetClusterBin()->GetUniqueCharge(this);
    }
    // Find the new weights for each link in this charge bin.
    for (ND::Share::TLinks::iterator clusterLink = GetClusterLinks().begin();
         clusterLink != GetClusterLinks().end(); ++clusterLink) {
        if (totalCharge>1E-9) {
            double w = (*clusterLink)->GetClusterBin()->GetUniqueCharge(this);
            (*clusterLink)->SetNewWeight(w/totalCharge);
        
        }
        else {
            (*clusterLink)->SetNewWeight(1.0/GetClusterLinks().size());
        }
    }
}

///////////////////////////////////////////////////////////////////////////
// TClusterBin
///////////////////////////////////////////////////////////////////////////

ND::Share::TClusterBin::TClusterBin() : fCluster(NULL) {
    fChargeLinks.clear();
}

ND::Share::TClusterBin::TClusterBin(ND::Share::TCluster* owner,
                                    ND::THandle<ND::TReconNode> node) 
    : fCluster(owner), fNode(node) {
    fChargeLinks.clear();
}

void ND::Share::TClusterBin::Dump() const {
    ND280Log("TClusterBin(" << std::hex << this << ")"
             << std::hex << " in TCluster(" << fCluster << ")"
             << std::dec << " w/ " << fChargeLinks.size() << " links");
#ifdef DUMP_LINKS
    ND::TND280Log::IncreaseIndentation();
    for (ND::Share::TLinks::const_iterator link = fChargeLinks.begin(); 
         link != fChargeLinks.end(); ++link) {
        (*link)->Dump();
    }
    ND::TND280Log::DecreaseIndentation();
#endif
}

double
ND::Share::TClusterBin::GetTotalCharge() const {
    double charge = 0;
    for (ND::Share::TLinks::const_iterator chargeLink 
             = GetChargeLinks().begin();
         chargeLink != GetChargeLinks().end(); ++chargeLink) {
        charge += (*chargeLink)->GetCharge();
    }
    return charge;
}

double
ND::Share::TClusterBin::GetUniqueCharge(
    const ND::Share::TChargeBin* cb) const {
    double charge = 0;
    for (ND::Share::TLinks::const_iterator chargeLink 
             = GetChargeLinks().begin();
         chargeLink != GetChargeLinks().end(); ++chargeLink) {
        if ((*chargeLink)->GetChargeBin() == cb) continue;
        charge += (*chargeLink)->GetCharge();
    }
    return charge;
}

/////////////////////////////////////////////////////////////////////
// TCluster
/////////////////////////////////////////////////////////////////////

ND::Share::TCluster::TCluster() {}

ND::Share::TCluster::TCluster(const ND::THandle<ND::TReconBase>& b) 
    : fReconObject(b) { }
    
ND::Share::TClusterBin& 
ND::Share::TCluster::AddClusterBin(ND::THandle<ND::TReconNode> node) {
    push_back(ND::Share::TClusterBin(this,node));
    return back();
}

void ND::Share::TCluster::Dump() const {
    ND280Log("TCluster(" << std::hex << this << ")"
             << " for " << ND::GetPointer(fReconObject)
             << std::dec << " w/ " << size() << " bins");
    ND::TND280Log::IncreaseIndentation();
    for (ND::Share::TCluster::const_iterator c = begin();
         c != end(); ++c) {
        c->Dump();
    }
    ND::TND280Log::DecreaseIndentation();
}

///////////////////////////////////////////////////////////////////
// TClusters
///////////////////////////////////////////////////////////////////
ND::Share::TClusters::TClusters() {}

ND::Share::TCluster& 
ND::Share::TClusters::AddCluster(const ND::THandle<ND::TReconBase>& b,
                                 ND::Share::TCharges& charges,
                                 ND::Share::TLinks& links) {
    push_back(ND::Share::TCluster(b));
    ND::Share::TCluster& cluster = back();

    for (ND::TReconNodeContainer::iterator n 
             = cluster.GetReconObject()->GetNodes().begin();
         n != cluster.GetReconObject()->GetNodes().end(); ++n) {
        cluster.AddClusterBin(*n);
    }

    // Create the links between the cluster bins and the charge bins.
    for (ND::Share::TCluster::iterator clusterBin = cluster.begin();
         clusterBin != cluster.end(); ++clusterBin) {
        ND::THandle<ND::THitSelection> hits 
            = clusterBin->GetNode()->GetObject()->GetHits();
        for (ND::THitSelection::iterator h = hits->begin();
             h != hits->end(); ++h) {
            ND::Share::TCharges::iterator c=charges.find((*h)->GetGeomId());
            if (c == charges.end()) {
                charges.insert(
                    ND::Share::TCharges::value_type(
                        (*h)->GetGeomId(), 
                        ND::Share::TChargeBin(*h)));
                c = charges.find((*h)->GetGeomId());
            }
            // Add a new link between the cluster and a charge bin.
            ND::Share::TLink* newLink 
                = new ND::Share::TLink(*clusterBin,c->second);
            links.push_back(newLink);
            // Correct for the attenuation.  This places all of the hits at
            // the center of the cluster (implicitly assuming that the
            // clusters are small in size).
            ND::THandle<ND::TShowerState> state 
                = clusterBin->GetNode()->GetState();
            double att = ND::TP0DAttenuationCurve::Get().GetAttenuation(
                (*h)->GetGeomId(),state->GetPosition().Vect());
            newLink->SetPhysicsWeight(1.0/att);
        }
    }

    return cluster;
}
