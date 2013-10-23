#ifndef TBestTube_hxx_seen
#define TBestTube_hxx_seen

#include <TReconCluster.hxx>

namespace CP {
    class TBestTube;
};

/// Take a set of clusters, and find a tube that contains a lot of them.  This
/// is kind of like a poor mans Hough transform, but has the advantage that it
/// works in 3 dimensions without exploding the memory usage.  It currently
/// uses a pure stochastic algorithm, but the Tabu search algorithm could also
/// be used. 
class CP::TBestTube {
public:
    /// Create a new tube search algorithm.
    TBestTube() {}
    ~TBestTube() {}

    /// Process a set of input clusters.  This takes an input
    /// TReconObjectContainer that is searched for the seed, and copies it to
    /// an internal area.  The FillSeed and FillRemains methods are then used
    /// to get the clusters that are part of the seed, and the clusters that
    /// were not used.
    void Process(const CP::TReconObjectContainer& input);

    /// Fill an output TReconObjectContainer with the clusters that are part
    /// of the seed.  This cannot be called before the Process method.
    void FillSeed(CP::TReconObjectContainer& seed);

    /// Fill an output TReconObjectContainer with the cluster that are NOT
    /// part of the seed.  This cannot be called before the Process method.
    void FillRemains(CP::TReconObjectContainer& remains);

private:
    /// Find a new iterator with a random step, but don't go outside of
    /// [begin,end).  This makes sure that the new iterator is valid.  The
    /// first argument, "i", is the current value of the iterator.  The extent
    /// of the container is specified by begin to end.
    CP::TReconObjectContainer::iterator 
    Randomize(CP::TReconObjectContainer::iterator i, 
              CP::TReconObjectContainer::iterator begin,
              CP::TReconObjectContainer::iterator end);

    /// Find the hits in a tube (as defined by it's end point positions).  The
    /// hits will be ordered from end1 to end2.
    void FindTube(const TVector3& end1, const TVector3& end2,
                  CP::TReconObjectContainer::iterator begin, 
                  CP::TReconObjectContainer::iterator end,
                  CP::TReconObjectContainer& output);

    /// Determine the "goodness" of a tube (as defined by it's end point
    /// positions).  A better tube will have a higher weight.
    double TubeWeight(const TVector3& end1, const TVector3& end2,
                      CP::TReconObjectContainer::iterator begin, 
                      CP::TReconObjectContainer::iterator end);

    /// Determine the "goodness" of a tube (as defined by the positions of the
    /// clusters at the end).  A better tube will have a higher weight.  This
    /// feeds the positions of the end clusters into the "other" version of
    /// TubeWeight.
    double TubeWeight(CP::TReconObjectContainer::iterator end1,
                      CP::TReconObjectContainer::iterator end2,
                      CP::TReconObjectContainer::iterator begin, 
                      CP::TReconObjectContainer::iterator end);

    /// The clusters that are copied so they can be returned via the GetTrack,
    /// FillRemains, and FillSeed methods.
    CP::TReconObjectContainer fClusterContainer;

    /// A flag that a valid tube has been found.
    bool fFoundTube;

    /// The position of the first end point of the tube.
    TVector3 fEnd1;

    /// The position of the second end point of the tube.
    TVector3 fEnd2;
};
#endif
