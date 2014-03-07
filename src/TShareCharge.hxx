#ifndef TP0DShareCharge_hxx_seen
#define TP0DShareCharge_hxx_seen

#include <iostream>
#include <list>

#include <TAlgorithm.hxx>
#include <TCaptLog.hxx>
#include <THandle.hxx>
#include <THit.hxx>

namespace CP {
    class TShareCharge;

    namespace ShareCharge {
        class TMeasurement;
        class TLink;
        class TMeasurementGroup;

        typedef std::list<TLink*> TLinks;
    };
    
};


/// An object describing a single measurement of the charge.  This measurement
/// may be shared by more than one TMeasurementGroup, and it's contribution is
/// going to be split over the TMeasurementGroups that contain it.  The
/// TMeasurement objects are "singles" in the since that the charge (or
/// energy) represented by a TMeasurement object is contained in exactly one
/// TMeasurement object (i.e. charge is not shared between TMeasurement
/// objects).  The TMeasurement object is connected to a single THit.
class CP::ShareCharge::TMeasurement {
public:
    typedef CP::THandle<CP::THit> Object;

    explicit TMeasurement(const Object& hit, double charge);

    /// Get the raw measurement object used to construct the TMeasurement.
    Object GetObject() const {return fObject;}

    /// Get the amount of charge for this measurement. 
    double GetCharge() const {return fCharge;}
    
    /// @{ Get the list of TMeasurementGroup objects that contain this
    /// TMeasurement object.
    TLinks& GetLinks() {return fLinks;}
    const TLinks& GetLinks() const {return fLinks;}
    /// @}

    /// Normalize the weights in the link list.
    void NormalizeWeights();

    /// Update the weights with the new weights in the link list.
    double UpdateWeights();

    /// Find the link weights for this measurement.
    void FindLinkWeights();

    /// Eliminate excess links.  This checks for any links that have almost
    /// zero weight and eliminates them.  When a link is eliminated, the
    /// measurement group that it links to have the charge set to zero.
    void EliminateLinks(double weightCut);

    /// Dump the measurement
    void Dump(bool dumpLinks = true) const;

private:
    /// The object associated with this measurement.  This is a single
    /// hit.
    Object fObject;

    /// The charge (i.e. the deposited energy) associated with this
    /// measurement.
    double fCharge;

    /// A list of links to the clusters bins which contain this measurement
    TLinks fLinks;
};

////////////////////////////////////////////////////////////////
/// An object describing a group of measurements that are combined into a
/// single cluster (or other "physics" related object).  This might represent
/// a track that could share measurements (i.e. hits) with another track, or
/// it could be a simple cluster.
///
/// In the context of CAPTAIN, this represents a 3D hit that is constructed
/// from 2D hits.  A single 2D hit (i.e. the measurement in the context of
/// CAPTAIN) might contribute to multiple 3D hits.
class CP::ShareCharge::TMeasurementGroup {
public:
    // The object referenced by the group.  In the context of CAPTAIN, this is
    // a 3D hit (i.e. a TReconHit).
    typedef CP::THandle<CP::THit> Object;

    /// Create a new measurement group and assign the owner.
    explicit TMeasurementGroup(CP::TShareCharge* owner, 
                               TMeasurementGroup::Object& object);

    /// Add a new measurement to this group.  This takes an object that will
    /// be linked to a TMeaasurement object (ie a THit), and the charge
    /// associated with that object.  If the input object is already
    /// associated with a TMeasurement, then a new link is added between that
    /// measurement and the current TMeasurementGroup.  If the input object is
    /// not associated with a TMeasurement, a new TMeasurement object is
    /// created and the link to this TMeasurementGroup is established.
    TMeasurement* AddMeasurement(TMeasurement::Object& object, double charge);

    /// @{ Get the list of measurements that are part of this cluster bin.
    const TLinks& GetLinks() const {return fLinks;}
    TLinks& GetLinks() {return fLinks;}
    // @}

    /// Get the node that is associated with this cluster bin.
    Object GetObject() const {return fObject;}

    /// Get the total charge in the group.  This returns the charge for the
    /// group adjusted by the current link weights (both physical and the
    /// weights being fitted by TShareCharge).
    double GetTotalCharge() const;

    /// Get the charge in the group not contributed by a particular
    /// measurement.
    double GetUniqueCharge(const TMeasurement* cb) const;

    void Dump(bool dumpLinks = true) const;

private:
    // This should never be used!
    TMeasurementGroup() :fOwner(NULL) {}

    /// The owner of this object;
    CP::TShareCharge* fOwner;

    /// The track or shower node associated with the TMeasurementGroup.
    Object fObject;

    /// The measurements that are part of this cluster bin.
    TLinks fLinks;
};

////////////////////////////////////////////////////////////////
/// A link between a TMeasurement object and one of the TMeasurementGroup
/// objects which contain it.  In addition to recording the connection between
/// a TMeasurement object and a TMeasurementGroup object, the link records the
/// fraction of the charge in the measurement that is associated with the
/// TMeasurementGroup.
class CP::ShareCharge::TLink {
public:
    TLink() : fWeight(1.0), fNewWeight(1.0), fPhysicsWeight(1.0), 
              fMeasurement(NULL), fMeasurementGroup(NULL) {}
        
    TLink(CP::ShareCharge::TMeasurementGroup* group, 
          CP::ShareCharge::TMeasurement* charge)
        : fWeight(1.0), fNewWeight(1.0), fPhysicsWeight(1.0), 
          fMeasurement(charge), fMeasurementGroup(group) { }

    /// Get the weight of the measurement in the linked measurement group.
    /// This is set as the result of the TShareCharge algoritm.
    double GetWeight() const {return fWeight;}

    /// Set the weight of the measurement in the linked measurement group.  
    void SetWeight(double w) {fWeight = w;}

    /// Get the new weight for the link.  This is used for internal
    /// bookkeeping during the calculation of the link weight by the
    /// TShareCharge algorithm.
    double GetNewWeight() const {return fNewWeight;}

    /// Set the new weight for the link.  This is used for internal
    /// bookkeeping during the calculation of the link weight by the
    /// TShareCharge algorithm.
    void SetNewWeight(double w) {fNewWeight = w;}

    /// Set the physics weight for this link.  The physics weight includes the
    /// effect of attenuation and other things that might need to be corrected
    /// when estimating how much a TMeasurement object contributes to a
    /// TMeasurementGroup.  This is not necessary in the context of CAPTAIN.
    void SetPhysicsWeight(double w) {fPhysicsWeight = w;}

    /// Get the physics weight for this link.
    double GetPhysicsWeight() const {return fPhysicsWeight;}

    /// Get the raw charge for the link.  This is the raw charge for the
    /// measurement without any weighting.
    double GetRawCharge() const {return GetMeasurement()->GetCharge();}

    /// Get the charge for the link.  This is the charge for the measurement
    /// corrected for any physics effects such as attenuation.
    double GetPhysicsCharge() const {return GetPhysicsWeight()*GetRawCharge();}

    /// Get the weighted charge that the measurement associated with this link
    /// adds to the measurement group associated with this link.
    double GetCharge() const {return GetWeight()*GetPhysicsCharge();}

    /// Get the measurement associated with the link.
    const CP::ShareCharge::TMeasurement* GetMeasurement() const {
        return fMeasurement;
    }

    /// Get the measurement group associated with this link.
    const CP::ShareCharge::TMeasurementGroup* GetGroup() const {
        return fMeasurementGroup;
    }

    /// Dump the values in the link.
    void Dump() const;

private:
    /// The weight of the link (between 0 and 1).  This gives the amount of
    /// charge in the measurement that should be added to the measurement group.
    double fWeight;

    /// The new weight of the link after relaxation.
    double fNewWeight;

    /// The physics based weighting between the measurement and the cluster
    /// bin.  This accounts for affects like attenuation and is a constant of
    /// the cluster to charge link.  The weight applied to the measurement
    /// should be actually fWeight*fPhysicsWeight.  The default value for the
    /// physics weight is 1.0.  Note that the physics weight can be greater
    /// than one.
    double fPhysicsWeight;

    /// The TMeasurement end of the link.
    CP::ShareCharge::TMeasurement* fMeasurement;

    /// The cluster end of the link;
    CP::ShareCharge::TMeasurementGroup* fMeasurementGroup;
};


/// Share charge among groups (ie TMeasurementGroup objects) of measurement
/// objects (ie TMeasurement objects) that might overlap.
///
/// The charge is shared among the groups by observing that if \f$n\f$
/// groups contribute \f$Q_1\f$, ... \f$Q_n\f$ to the total charge in the
/// event, then the total charge in a group \f$i\f$ can be expressed as
///
/// \f[ Q_i = w_i q_s + Q'_i \f] 
///
/// where \f$q_s\f$ is an element of charge shared between multiple groups,
/// \f$w_i\f$ is the fraction of the shared charge contributing to group
/// \f$i\f$, and \f$Q'_i\f$ is the charge that is not shared with other
/// groups.  For groups, the measurement element is taken to be a single hit
/// and it will typically be shared between one (i.e. not shared), two or
/// three groups.  For clarity, the following assumes the hit is shared
/// between three groups, but the derived result is general.
///
/// When three groups overlap and share a measurement \f$q_s\f$, the
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
/// share measurements with other clusters, the solution for all of the
/// weights contributing to an event must be solved iteratively.
///
/// This code is derived from the IMB3 event reconstruction.
///
/// \code
/// /* #Id: share_flux.c,v 1.1 1994/02/28 19:33:15 clark Exp mcgrew # */
/// /* Find the fraction of the total Cherenkov energy that is in each track. */
/// /* This is called by a routine that finds the flux as a function of */
/// /* direction.  The main entry point is (share_flux).   */
/// \endcode
class CP::TShareCharge {
public:
    friend class CP::ShareCharge::TMeasurementGroup;

    typedef std::list<CP::ShareCharge::TLink> Links;
    typedef std::list<CP::ShareCharge::TMeasurement> Measurements;
    typedef std::list<CP::ShareCharge::TMeasurementGroup> Groups;

    TShareCharge();
    virtual ~TShareCharge();

    /// Add a new measurement group.  This allocates the new group empty and
    /// then returns a reference.  The group will then need to have the
    /// measurements that make up the group added to it.
    CP::ShareCharge::TMeasurementGroup& 
    AddGroup(CP::ShareCharge::TMeasurementGroup::Object& object);

    void DumpGroups(bool dumpLinks = true) const;
    void DumpMeasurements(bool dumpLinks = true) const;
    void Dump(bool dumpLinks = true) const;

    /// Solve the coupled equations to find the optimal set of weight to share
    /// the charge measurements among the measurement groups.  After this has
    /// been called, the charge in the measurement groups have been updated.
    /// The couple equations are solved using interative relaxation.  The
    /// return value is the change in the last iteration. 
    double Solve(double tolerance = 1E-3, int iterations = 5000);

    /// Return the measurement groups.  This is how the result of the charge
    /// sharing is accessed.
    const Groups& GetGroups() const {return fGroups;}

    /// Set the weight cut used to remove links where the weight is
    /// effectively zero.
    void SetWeightCut(double w) {fWeightCut = w;}

    /// Get the weight cut.
    double GetWeightCut() const {return fWeightCut;}

private:
    /// This returns how much the weights have changed during the iteration.
    /// It should be called until the change is small.
    double RelaxWeights();

    /// Get an measurement from the collection of measurements.  If the
    /// measurement does not exist, then it will be added to the collection.
    /// This returns a reference to the measurement.  This is used by the
    /// TMeasurementGroup class (which is a friend).  The TMeasurementGroup
    /// object will need to add the appropriate links between itself and the
    /// TMeasurement object.
    CP::ShareCharge::TMeasurement* FindMeasurement(
        CP::ShareCharge::TMeasurement::Object& object, double charge);


    /// Create a new link between a measurement and a measurement group.  THis
    /// is used by the TMeasurementGroup object.
    CP::ShareCharge::TLink* CreateLink(
        CP::ShareCharge::TMeasurementGroup* group,
        CP::ShareCharge::TMeasurement* measurement);
    
    /// All of the measurements associated with this object.
    Measurements fMeasurements;

    /// All of the measurement groups associated with this object.
    Groups fGroups;

    /// All of the links between measurements and groups associated with this
    /// object.
    Links fLinks;

    /// The charge fraction below which links are eliminated.  This is used to
    /// remove links where the weight has effectively gone to zero.
    double fWeightCut;

};



#endif
