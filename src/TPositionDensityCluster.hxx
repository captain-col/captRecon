#ifndef TPositionDensityCluster_seen
#define TPositionDensityCluster_seen

#include "TIterativeNeighbors.hxx"

#include "TCaptLog.hxx"

#include <vector>
#include <list>
#include <functional>
#include <algorithm>

namespace CP {
    template <class PositionHandle> class TPositionDensityCluster;
}

/// A class that performs density-based clustering using the DBSCAN algorithm
/// using a cartesian metric.  It works on any class that implements the
/// PositionHandle concept.  Here is a code snippet using the this template:
///
/// \code
///    // Make a typedef for the ClusterAlgorithm.
///    const double maxDist = 20*unit::cm;
///    const unsigned int minPoints = 4;
///    
///    std::auto_ptr<TPositionDensityCluster<CP::THandle<CP::THit> >
///        xcluster(new TPositionDensityCluster(minPoints, maxDist));
///    xcluster->Cluster(xHits);
/// \endcode
///
/// For a detailed description of the density-based clustering, Google
/// keyword: density-based clustering.  This has an interface similar to the
/// TTmplDensityCluster template which can handle generalized distance
/// definitions.  The TPositionDensityCluster class is limited to clustering
/// using the cartesian distances between hit positions.
///
/// The template argument must implement the PositionHandle concept which
/// means it behaves as a pointer to a class that provides the GetPosition()
/// method which typically returns a TVector3, or TLorentzVector.  The
/// requirement is that the objectd returned by the GetPosition() method must
/// be returned by reference, and must implement the X(), Y() and Z() methods
/// returning a float or double value.
template <class PositionHandle>
class CP::TPositionDensityCluster {
public:
    /// A collection of points for use in clustering.  A Points collection is
    /// returned by GetCluster().
    typedef std::list<PositionHandle> Points;

    /// A type to be used with the neighbor position clustering.  The first
    /// element is a "color" for the entry, and the second is the handle that
    /// gets returned to the caller.
    struct NeighborEntry {
        NeighborEntry(): fColorIndex(-1) {}
        NeighborEntry(int c,const PositionHandle& h) 
            : fColorIndex(c), fHandle(h) {}
        int fColorIndex;
        PositionHandle fHandle;
    };

    /// The k-d tree neighbors.
    typedef typename CP::TIterativeNeighbors<NeighborEntry> Neighbors;

    /// Create a density clustering class that requires at least minPts within
    /// a distance of maxDist.  The distance is the cartesion distance between
    /// the points.
    TPositionDensityCluster(std::size_t minPts, double maxDist);
    virtual ~TPositionDensityCluster() {}
    
    /// Cluster a group of objects between the begin and end iterator.  The
    /// results are accessed using GetCluster().  The iterators must implement
    /// the concept of an iterator to objects of type PositionHandle.
    template <class InputIterator>
    void Cluster(InputIterator begin, InputIterator end);

    /// Return the number of clusters found by the density clustering.  This
    /// is only valid after the Cluster() method has been used.
    std::size_t GetClusterCount() { return fClusters.size(); }

    /// Get the i-th cluster.  This is only valid after the cluster method has
    /// been used.  If the index is equal to the number of found clusters,
    /// then the return value will be the list of unclustered points.
    const Points& GetCluster(std::size_t i) const {
        if (fClusters.size() <= i) return fRemaining;
        return fClusters.at(i); 
    }
    
protected:

    /// Find the set of points with highest density in input.  If the density
    /// is greater that fMinPoints, then return the set of points in output.
    /// Points that are in the seed are added to the output using GetNeighbors
    /// and have their color changed from black to white.  This returns true
    /// if a new cluster seed was found.
    bool FindSeeds(Neighbors& input, Points& output);

    /// Find the neighbors for a reference hit in the input points and place
    /// them into the output.  Neighbors are defined as all points for which
    /// the distance to the reference is less than fMaxDist (a point is a
    /// neighbor to itself).  Any neighbors that are not currently in a
    /// cluster will be returned in the output variable.  Any point added to
    /// the output has it's color changed from black to white.
    void GetNeighbors(PositionHandle reference, 
                      Neighbors& input, Points& output);

    /// Count the neighbors for a reference in the input points, but do not
    /// return a copy of the neighbors.  The point is not counted as a
    /// neighbor.  This counts points even if they have already been added to
    /// a cluster.
    std::size_t CountNeighbors(PositionHandle reference, Neighbors& input);

    /// Color all of the grey handles in the map to black.
    void MakeGreyBlack();

    /// Color all of the grey handles in the map to white.
    void MakeGreyWhite();

private:
    /// The minimum number of points that must be within the fMaxDist
    /// radius of the current point.  If there are at least fMinPoints with in
    /// the fMaxDist of the current point, the cluster will be expanded by
    /// the neighbors.  The distance between points is defined by
    /// fMetricModel.
    unsigned int fMinPoints;

    /// The maximum distance between points for which points are defined as
    /// being neighbors. 
    double fMaxDist;

    /// The clusters that have been found.
    std::vector<Points> fClusters;

    /// The objects that didn't make it into a cluster.
    Points fRemaining;

    /// The mapping of color to each point in the neighborTree.  This is
    /// indexed by the fColorIndex field of the NeighborEntry.  The color map
    /// has three shades: 0) white -- this is when the point is free.  1) grey
    /// -- this is when a point is currently being checked.  2) black -- this
    /// is when a point has been included in a cluster.
    enum {kWhite=0, kGrey, kBlack};
    std::vector<int> fColorMap;

    /// The sorting comparison to order the clusters at the end of the search.
    template <class T>
    struct ClusterOrdering :
        public std::binary_function <T, T, bool > {
        bool operator() (const T& lhs, const T& rhs) {
            return lhs.size() > rhs.size();
        }
    };
};

////////////////////////////////////////////////////////////////
// Define the TPositionDensityCluster class methods.
////////////////////////////////////////////////////////////////

template <class PositionHandle>
CP::TPositionDensityCluster<PositionHandle>::TPositionDensityCluster(
    std::size_t MinPts, double maxDist) 
    : fMinPoints(MinPts), fMaxDist(maxDist) { }

template <class PositionHandle>
void CP::TPositionDensityCluster<PositionHandle>::MakeGreyBlack() {
    for (std::vector<int>::iterator c = fColorMap.begin();
         c != fColorMap.end(); ++c) {
        if (*c == kGrey) *c = kBlack;
    }
}

template <class PositionHandle>
void CP::TPositionDensityCluster<PositionHandle>::MakeGreyWhite() {
    for (std::vector<int>::iterator c = fColorMap.begin();
         c != fColorMap.end(); ++c) {
        if (*c == kGrey) *c = kWhite;
    }
}

template <class PositionHandle>
template <class InputIterator>
void CP::TPositionDensityCluster<PositionHandle>::Cluster(
    InputIterator begin, InputIterator end) {
    
    // Clear out the  internal data structures.
    fClusters.clear();

    // Insert the input into the tree.
    Neighbors neighborTree;
    int index = 0;
    for (InputIterator handle = begin; handle != end; ++handle) {
        neighborTree.AddPoint(NeighborEntry(index++,*handle),
                              (*handle)->GetPosition().X(),
                              (*handle)->GetPosition().Y(),
                              (*handle)->GetPosition().Z());
    }
    fColorMap.resize(index);
    std::fill(fColorMap.begin(), fColorMap.end(), kWhite);

    CaptNamedDebug("cluster", "Input points: " << index);

    // Now continue removing points until there aren't any more points, or a
    // seed isn't found.
    while (true) {
        Points seeds; 

        // Find the next set of seeds to start a cluster.  The seeds are
        // enough to form a cluster, and are marked in FindSeeds (by way of
        // GetNeighbors) as being in a cluster.
        seeds.clear();
        if (!FindSeeds(neighborTree,seeds)) {
            CaptNamedDebug("cluster", "No seed found");
            break;
        }

        CaptNamedDebug("cluster", "Start seed with " 
                       << seeds.size() << " points");

        // Start a new cluster and add the current seeds to it.
        Points cluster;
        std::copy(seeds.begin(), seeds.end(), std::back_inserter(cluster));

        Points tmp;
        while (!seeds.empty()) {
            PositionHandle handle = seeds.front();
            seeds.erase(seeds.begin());
            std::size_t count = CountNeighbors(handle, neighborTree);
            if (count < fMinPoints) continue;
            tmp.clear();
            GetNeighbors(handle,neighborTree,tmp);
            MakeGreyBlack();
            std::copy(tmp.begin(),tmp.end(),std::back_inserter(seeds));
            std::copy(tmp.begin(),tmp.end(),std::back_inserter(cluster));
        }

        CaptNamedDebug("cluster", "Cluster with "
                       << cluster.size() << " points");
        fClusters.push_back(cluster);
    }

    std::sort(fClusters.begin(), fClusters.end(),
              ClusterOrdering<Points>());

    fRemaining.clear();
    for (typename Neighbors::value_iterator p = neighborTree.begin_values(); 
         p != neighborTree.end_values(); ++p) {
        if (fColorMap[p->fColorIndex]) continue;
        fRemaining.push_back(p->fHandle);
    }
    
}
    
template <class PositionHandle>
bool CP::TPositionDensityCluster<PositionHandle>::FindSeeds(
    Neighbors& in, Points& out) {
    out.clear();
    double dist2 = fMaxDist*fMaxDist;
    int seedCount = 0;
    for (typename Neighbors::value_iterator p = in.begin_values(); 
         p != in.end_values(); ++p) {
        if (fColorMap[p->fColorIndex] == kBlack) {
            continue;
        }
        std::size_t count = 0;
        PositionHandle h = p->fHandle;
        typename Neighbors::iterator neighbor = in.begin(h->GetPosition().X(),
                                                h->GetPosition().Y(),
                                                h->GetPosition().Z());
        typename Neighbors::iterator neighbor_end = in.end();
        for (;neighbor != neighbor_end; ++neighbor) {
            // Stop looking at distant handles.
            if (neighbor->second > dist2) {
                break;
            }
            // Don't count a handle that's already in a cluster.
            if (fColorMap[neighbor->first.fColorIndex] == kBlack) {
                continue;
            }
            ++count;
            // Short-circuit the search.
            if (count>2*fMinPoints) break;
        }
        // Increment count since the current point is also part of the seed.
        ++count;
        // Not enough points, so look at the next value.
        if (count<fMinPoints) {
            continue;
        }
        // We already have a bigger seed.
        if (count<out.size()) {
            continue;
        }
        // We found that "h" is a better seed, so copy it into out. This
        // adds the test hit as well as it's neighbors to the seeds.
        out.clear();
        MakeGreyWhite();
        GetNeighbors(h, in, out);
        // Sort-circuit the search if we've found a good enough seed.
        if (seedCount > 5 && count > 3*fMinPoints) break;
        ++seedCount;
    }
    if (out.size() < fMinPoints) {
        MakeGreyWhite();
        return false;
    }
    MakeGreyBlack();
    return true;
}

template <class PositionHandle>
void CP::TPositionDensityCluster<PositionHandle>::GetNeighbors(
    PositionHandle pnt, Neighbors& in, Points& out) {
    double dist2 = fMaxDist*fMaxDist;
    typename Neighbors::iterator neighbor = in.begin(pnt->GetPosition().X(),
                                                     pnt->GetPosition().Y(),
                                                     pnt->GetPosition().Z());
    typename Neighbors::iterator neighbor_end = in.end();
    for (;neighbor != neighbor_end; ++neighbor) {
        // Stop looking at distant handles.
        if (neighbor->second > dist2) break;
        // Don't add handle that's already in a cluster.
        if (fColorMap[neighbor->first.fColorIndex]) continue;
        // Mark the handle as in a cluster.
        fColorMap[neighbor->first.fColorIndex] = kGrey;
        // Add the value to the output
        out.push_back(neighbor->first.fHandle);
    }
}
    
template <class PositionHandle>
std::size_t CP::TPositionDensityCluster<PositionHandle>::CountNeighbors(
    PositionHandle pnt, Neighbors& in) {
    double dist2 = fMaxDist*fMaxDist;
    std::size_t count = 0;
    typename Neighbors::iterator neighbor = in.begin(pnt->GetPosition().X(),
                                                     pnt->GetPosition().Y(),
                                                     pnt->GetPosition().Z());
    typename Neighbors::iterator neighbor_end = in.end();
    for (;neighbor != neighbor_end; ++neighbor) {
        // Stop looking at distant handles.
        if (neighbor->second > dist2) break;
        // Don't count the current handle.
        if (neighbor->second < 1E-6) continue;
        ++count;
    }
    return count;
}
#endif

