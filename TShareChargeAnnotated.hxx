#ifndef TP0DShareCharge_hxx_seen
#define TP0DShareCharge_hxx_seen

#include <iostream>
#include <set>
#include <list>
#include <map>

#include <TAlgorithm.hxx>
#include <TND280Log.hxx>
#include <TReconBase.hxx>
#include <TReconNode.hxx>
#include <TAlgorithmResult.hxx>
#include <TGeometryId.hxx>

namespace ND {
    namespace Share {
        class TChargeBin;
        class TLink;
        class TClusterBin;
        class TCluster;
        class TClusters;

        /// A list of links between clusters and charge bins.
        typedef std::list<TLink*> TLinks;

        /// A set of all TChargeBins.
        typedef std::map<ND::TGeometryId, ND::Share::TChargeBin> TCharges;
    };
    
    class TP0DShareCharge;
};

////////////////////////////////////////////////////////////////
/// An object describing the charge in a single bin that is part of a cluster
/// of hits.  This collects the charge from a single hit.  The charge bin will
/// often (usually?) contribute to more than one TClusterBin.
class ND::Share::TChargeBin {
public:
    explicit TChargeBin(const ND::THandle<ND::THit>& hit);

    /// Get the amount of charge in this bin.
    double GetCharge() const {return fHit->GetCharge();}
    
    /// @{ Get the list of cluster bins that contains this charge bin.
    TLinks& GetClusterLinks() {return fClusterLinks;}
    const TLinks& GetClusterLinks() const {return fClusterLinks;}
    /// @}

    /// Normalize the weights in the link list.
    void NormalizeWeights();

    /// Find the link weights for this charge bin.
    void FindLinkWeights();

    /// Update the weights with the new weights in the link list.
    double UpdateWeights();

    void Dump() const;

private:
    /// The hits associated with this charge bin.  This is normally a single
    /// hit.
    ND::THandle<ND::THit> fHit;

    /// A list of links to the clusters bins which contain this charge bin
    TLinks fClusterLinks;
};

////////////////////////////////////////////////////////////////
/// An object describing a cluster of hits.  This is associated with a node in
/// the parent reconstruction object.  This usually represents a single
/// TReconCluster, but could represent any collection of hits.  
class ND::Share::TClusterBin {
public:
    TClusterBin();

    /// Create a new cluster bin and assign who owns it.
    explicit TClusterBin(ND::Share::TCluster* owner, 
                         ND::THandle<ND::TReconNode> node);

    /// Get the cluster that contains this bin.
    const TCluster* GetCluster() const {return fCluster;}

    /// @{ Get the list of charge bins (e.g. layers or hits) that are part of
    /// this cluster bin.
    const TLinks& GetChargeLinks() const {return fChargeLinks;}
    TLinks& GetChargeLinks() {return fChargeLinks;}
    // @}

    /// Get the node that is associated with this cluster bin.
    ND::THandle<ND::TReconNode> GetNode() const {return fNode;}

    /// Get the total charge in the cluster.  This returns the charge for the
    /// cluster adjusted by the current link weights (both physical and the
    /// weights being fitted by TP0DShareCharge).
    double GetTotalCharge() const;

    /// Get the charge in the cluster not contributed by a particular charge
    /// bin.
    double GetUniqueCharge(const TChargeBin* cb) const;

    void Dump() const;

private:

    /// The cluster that contains this bin.
    TCluster* fCluster;

    /// The track or shower node associated with the TClusterBin.
    ND::THandle<ND::TReconNode> fNode;

    /// The charge bins that are part of this cluster bin.
    TLinks fChargeLinks;
};

////////////////////////////////////////////////////////////////
/// An object describing a cluster that is broken into cluster bins.  This is
/// really just a shower or track object broken up into it's nodes
/// (TClusterBins).
class ND::Share::TCluster: public std::list<ND::Share::TClusterBin> {
public:
    TCluster();
    explicit TCluster(const ND::THandle<ND::TReconBase>& b);
    
    /// Add a new cluster bin to the cluster which is then returned so that
    /// new charge bins can be added to the cluster bin.  This is the generic
    /// low level way to add new cluster bins to a TCluster.
    ND::Share::TClusterBin& AddClusterBin(ND::THandle<ND::TReconNode> node);

    /// Get the reconstruction object that is associated with this cluster.
    ND::THandle<ND::TReconBase> GetReconObject() const {return fReconObject;}

    void Dump() const;

private:
    /// The reconstruction object that is represented by this cluster.
    ND::THandle<ND::TReconBase> fReconObject;
};

class ND::Share::TClusters: public std::list<ND::Share::TCluster> {
public:
    TClusters();
    
    /// Add a new cluster to the vector of clusters.
    ND::Share::TCluster& AddCluster(const ND::THandle<ND::TReconBase>& b,
                                    ND::Share::TCharges& charges,
                                    ND::Share::TLinks& links);
};

////////////////////////////////////////////////////////////////
/// A link between a TChargeBin object and a TClusterBin.
class ND::Share::TLink {
public:
    TLink() : fWeight(1), fNewWeight(1), fPhysicsWeight(1), 
              fChargeBin(NULL), fClusterBin(NULL) {}
        
    TLink(ND::Share::TClusterBin& cluster, ND::Share::TChargeBin& charge)
        : fWeight(1), fNewWeight(1), fPhysicsWeight(1), 
          fChargeBin(&charge), fClusterBin(&cluster) {
        charge.GetClusterLinks().push_back(this);
        cluster.GetChargeLinks().push_back(this);
    }

    /// Get the weight of the charge bin in the linked cluster bin.
    double GetWeight() const {return fWeight;}

    /// Set the weight of the charge bin in the linked cluster bin.
    void SetWeight(double w) {fWeight = w;}

    /// Get the new weight for the link.
    double GetNewWeight() const {return fNewWeight;}

    /// Set the new weight for the link.
    void SetNewWeight(double w) {fNewWeight = w;}

    /// Set the physics weight for this link.  The physics weight includes the
    /// effect of attenuation and other things that need to be corrected when
    /// estimating how much a TChargeBin object contributes to a TClusterBin.
    void SetPhysicsWeight(double w) {fPhysicsWeight = w;}

    /// Get the physics weight for this link.
    double GetPhysicsWeight() const {return fPhysicsWeight;}

    /// Get the raw charge for the link.  This is the raw charge for the
    /// charge bin without any weighting.
    double GetRawCharge() const {return GetChargeBin()->GetCharge();}

    /// Get the charge for the link.  This is the charge for the charge bin
    /// corrected for any physics effects such as attenuation.
    double GetPhysicsCharge() const {return GetPhysicsWeight()*GetRawCharge();}

    /// Get the weighted charge for the link.
    double GetCharge() const {return GetWeight()*GetPhysicsCharge();}

    /// Get the charge bin associated with the link.
    const ND::Share::TChargeBin* GetChargeBin() const {return fChargeBin;}

    /// Get the cluster bin associated with this link.
    const ND::Share::TClusterBin* GetClusterBin() const {return fClusterBin;}

    /// Dump the values in the link.
    void Dump() const {
        ND280Log("TLink(" << std::hex << this << ")"
                 << std::dec <<std::setprecision(3) << " w " << fWeight
                 << std::dec <<std::setprecision(3) << " n " << fNewWeight
                 << std::dec <<std::setprecision(3) << " p " << fNewWeight
                 << std::hex << " c " << fClusterBin
                 << std::hex << " q " << fChargeBin
                 << std::dec);
    }

private:
    /// The weight of the link (between 0 and 1).  This gives the amount of
    /// charge in the charge bin that should be added to the cluster bin.
    double fWeight;

    /// The new weight of the link after relaxation.
    double fNewWeight;

    /// The physics based weighting between the charge bin and the cluster
    /// bin.  This accounts for affects like attenuation and is a constant of
    /// the cluster to charge link.  The weight applied to the charge bin
    /// should be actually fWeight*fPhysicsWeight.  The default value for the
    /// physics weight is 1.0.  Note that the physics weight can be greater
    /// than one.
    double fPhysicsWeight;

    /// The charge end of the link.
    ND::Share::TChargeBin* fChargeBin;

    /// The cluster end of the link;
    ND::Share::TClusterBin* fClusterBin;
};


/// Share charge among TReconShower objects that might overlap.  The input
/// algorithm result should contain a recon object container with vertices
/// containing showers and TReconPID objects.  The TReconPID objects will
/// have been derived from the TP0DPairwiseVertex_PID algorithm and are
/// ignored by this algorithm.  If there are no showers, then the input
/// vertex is just copied to the output.
///
/// The charge is shared among the showers by observing that if \f$n\f$
/// showers contribute \f$Q_1\f$, ... \f$Q_n\f$ to the total charge in the
/// event, then the total charge in a shower \f$i\f$ can be expressed as
///
/// \f[ Q_i = w_i q_s + Q'_i \f] 
///
/// where \f$q_s\f$ is an element of charge shared between multiple showers,
/// \f$w_i\f$ is the fraction of that charge generated by the shower, and
/// \f$Q'_i\f$ is the charge that is not in the shared element.  For showers,
/// the shared charge element is taken to be a single hit and it will
/// typically be shared between one (i.e. not shared), two or three showers.
/// For clarity, the following assumes the hit is shared between three
/// showers, but the derived result is general.
///
/// When three showers overlap and share an element of charge \f$q_s\f$, the
/// ratio of the weights are related.
///
/// \f[\frac{w_1}{w_2} = \frac{Q'_1}{Q'_2},\ \frac{w_1}{w_3} = \frac{Q'_1}{Q'_3},\ \frac{w_2}{w_3} = \frac{Q'_2}{Q'_3}\f]
///
/// or equivalently
///
/// \f[\frac{w_1}{Q'_1} = \frac{w_2}{Q'_2},\ \frac{w_1}{Q'_1} = \frac{w_3}{Q'_3},\ \frac{w_2}{Q'_2} = \frac{w_3}{Q'_3}\f]
///
/// subject to the constraint that \f$w_1 + w_2 + w_3 = 1\f$.
///
/// The optimal choice of weights can then be found by minimizing
///
/// \f[S = \left(\frac{w_1}{Q'_1}-\frac{w_2}{Q'_2}\right)^2 + \left(\frac{w_1}{Q'_1} - \frac{w_3}{Q'_3}\right)^2 + \left(\frac{w_2}{Q'_2} - \frac{w_3}{Q'_3}\right)^2 + \left(w_1 + w_2 + w_3 - 1 \right)^2\f]
///
/// with respect to \f$w_1\f$, \f$w_2\f$, and \f$w_3\f$.  This has the simple
/// solution that
///
/// \f[w_1 = \frac{Q'_1}{Q'_1 + Q'_2 + Q'_3}\f]
///
/// with similar solutions for all weights.
///
/// Finally, since the \f$Q'_i\f$ entering into the solution for \f$w_i\f$ may
/// share other charge elements the solution for all of the weights
/// contributing to an event must be solved iteratively.
///
/// This code is derived from the IMB3 event reconstruction.
///
/// \code
/// /* #Id: share_flux.c,v 1.1 1994/02/28 19:33:15 clark Exp mcgrew # */
/// /* Find the fraction of the total Cherenkov energy that is in each track. */
/// /* This is called by a routine that finds the flux as a function of */
/// /* direction.  The main entry point is (share_flux).   */
/// \endcode
class ND::TP0DShareCharge: public ND::TAlgorithm {
public:
    TP0DShareCharge();
    virtual ~TP0DShareCharge();
    
    ND::THandle<ND::TAlgorithmResult> Process(const ND::TAlgorithmResult&);

private:
    /// Actually share the charge between showers.  The showers should
    /// contained in a TReconObjectContainer.  A new container is created with
    /// new showers.
    ND::THandle<ND::TReconObjectContainer> 
    ShareCharge(const ND::TReconObjectContainer& input) const;

    /// Do one step of relaxation for all of the link weights between
    /// overlapping cluster bins and charge bins.
    double RelaxWeights(ND::Share::TCharges& charges) const;
};
#endif
