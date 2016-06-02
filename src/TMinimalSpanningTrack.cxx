#include "TMinimalSpanningTrack.hxx"
#include "ClusterDistance.hxx"
#include "CreateTrack.hxx"

#include <THandle.hxx>
#include <TReconTrack.hxx>
#include <TReconCluster.hxx>
#include <TCaptLog.hxx>
#include <HEPUnits.hxx>
#include <TUnitsTable.hxx>
#include <TRuntimeParameters.hxx>

#include <TRandom.h>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/betweenness_centrality.hpp>

#include <memory>
#include <cmath>

// A namespace for the type definitions needed to work with the Boost Graph
// Library.  These could be made into local typedefs, or they could be put
// into the TMinimalSpanningTrack class, but they are only used in this file.
// Notice that "vertex" in this code refers to the graph theory "vertex", and
// not the event "vertex".
namespace MST {
    /// This is a structure of any extra fields that need to be attached to
    /// the vertex.
    struct vertex_properties_t {
        /// The is cluster associated with this vertex.
        CP::THandle<CP::TReconCluster> cluster;
        /// This is the parent vertex of this one.  There is no parent if this
        /// is less than zero.
        int parent;
        /// This will be filled with the distance to the parent.
        double distance;
        /// The total distance to the root.
        double rootDistance;
        /// This is the number of steps from the root.
        int steps;
        /// This is the number of children for this vertex.
        int children;
        /// This is set to true after the vertex has it's cluster added to a
        /// track. 
        bool marked;
    };
    
    /// This is a structure of any extra fields that need to be attached to
    /// the edge.  This example is using Dijkstra's algorithm to find a MST,
    /// so the edge needs to have a length (or weight).  That parameter is
    /// declared.  Any extra fields could be added.  This length will be
    /// declared to Dijkstras algorithm using named parameters
    /// (i.e. boost::weight_map()) and boost::get()). This works out to be
    /// boost::weight_map(boost::get(&edge_properties_t::length,g));
    struct edge_properties_t {
        /// The length of the edge.
        double length;
    };
    
    /// The graph type being used in this analysis.  This says that the edges
    /// and vertices are going to be stored in a std::vector, the graph is
    /// undirected, and that the vertices and edges have extra fields defined
    /// in the vertex_properties_t and edge_properties_t structures.
    typedef boost::adjacency_list < 
        boost::vecS,              /// The edge container type
        boost::vecS,              /// The vertex container type.
        boost::undirectedS,       /// The graph is undirected.
        MST::vertex_properties_t, /// The extra fields for the vertices.
        MST::edge_properties_t    /// The extra fields for the edges.
        > graph_t;
    
    /// This is the "index" into the graph for the vertices.  It's used as 
    ///
    /// \code
    /// graph_t g;
    /// vertex_t v;
    /// g[v].distance = blah; // access the vertex distance field.
    /// \endcode
    typedef boost::graph_traits < graph_t >::vertex_descriptor vertex_t;

    /// An iterator over the vertices in the graph.
    typedef boost::graph_traits< graph_t >::vertex_iterator vertex_iterator_t;

    /// The "index" into the graph for the edges.
    typedef boost::graph_traits < graph_t >::edge_descriptor edge_t;

    typedef boost::graph_traits< graph_t >::adjacency_iterator adjacency_iterator_t;
};

CP::TMinimalSpanningTrack::TMinimalSpanningTrack()
    : TAlgorithm("TMinimalSpanningTrack", 
                 "Find Tracks Using on a Minimal Spanning Tree") {
}

CP::TMinimalSpanningTrack::~TMinimalSpanningTrack() { }

CP::THandle<CP::TAlgorithmResult>
CP::TMinimalSpanningTrack::Process(const CP::TAlgorithmResult& input,
                           const CP::TAlgorithmResult&,
                           const CP::TAlgorithmResult&) {

    CP::THandle<CP::TReconObjectContainer> inputObjects 
        = input.GetResultsContainer();

    CaptLog("TMinimalSpanningTrack Process " << GetEvent().GetContext());

    if (!inputObjects) {
        CaptError("No input objects");
        return CP::THandle<CP::TAlgorithmResult>();
    }

    // Create the output containers.
    CP::THandle<CP::TAlgorithmResult> result = CreateResult();
    std::unique_ptr<CP::TReconObjectContainer> 
        final(new CP::TReconObjectContainer("final"));

    // Find the clusters in the input object and copy them into the
    // remainingClusters object.  Any non-cluster objects are copied directly
    // to final.
    CP::TReconObjectContainer remainingClusters;
    for (CP::TReconObjectContainer::iterator i = inputObjects->begin();
         i != inputObjects->end(); ++i) {
        CP::THandle<CP::TReconCluster> cluster = *i;
        if (!cluster) {
            final->push_back(*i);
            continue;
        }
        remainingClusters.push_back(*i);
    }

    // Loop until all of the remaining clusters are handled.
    int throttle = remainingClusters.size();
    do {
        CaptNamedLog("MST","Remaining iterations " << throttle 
                  << " Remaining Clusters " << remainingClusters.size());

        // Create the new graph.
        MST::graph_t g(remainingClusters.size());

        // Insert the remaining clusters into the graph.  A new graph is
        // created for each iteration.
        MST::vertex_iterator_t vi, vi_end;
        for (boost::tie(vi,vi_end) = boost::vertices(g); 
             vi != vi_end; ++vi) {
            std::size_t index = *vi;
            g[*vi].cluster = remainingClusters[index];
            g[*vi].parent = -1;
            g[*vi].distance = 0;
            g[*vi].rootDistance = 0;
            g[*vi].steps = 0;
            g[*vi].children = 0;
            g[*vi].marked = false;
        }

        // Add the edges.  The distance cut is a preliminary cut to reduce the
        // total number of edges.  The hit distance cut is the actual closest
        // approach between the hits in the clusters.  These should be made
        // into parameters.
        double fDistCut = 300*unit::mm;
        double fHitDistCut = 100*unit::mm;
        for (boost::tie(vi,vi_end) = boost::vertices(g); vi != vi_end; ++vi) {
            MST::vertex_iterator_t vj, vj_end;
            for (vj = vi+1; vj != vi_end; ++vj) {
                double dist = (g[*vi].cluster->GetPosition().Vect()
                               - g[*vj].cluster->GetPosition().Vect()).Mag();
                if (dist > fDistCut) continue;
                double hitDist = CP::ClusterVicinity(*g[*vi].cluster,
                                                     *g[*vj].cluster);
                if (hitDist > fHitDistCut) continue;
                hitDist = CP::ClusterDistance(*g[*vi].cluster,
                                              *g[*vj].cluster);
                MST::edge_t edge = boost::add_edge(*vi, *vj, g).first;
                g[edge].length = hitDist;
            }
        }

        // There are not enough edges!  Stop now.
        if (boost::num_edges(g) < 2) break;

        int bestRoot = -1;
#ifdef USE_MOST_CENTRAL
        // Find the most central vertex on the principle that it will be close
        // to the event vertex.  This is going to be the root.  The Brandes
        // betweenness centrality is O(V*E) [V is number of vertices, and E is
        // number of edges], so don't run it when there are too many edges or
        // vertices.  When there are to many vertices, use a randomly chosen
        // vertex (assuming it will be relatively central to the position
        // distribution).  The actual choice of "most central" is not
        // critical.
        std::size_t mostCentral = gRandom->Integer(boost::num_vertices(g));
        double maxCentrality = -1;
        if (boost::num_vertices(g)*boost::num_edges(g) < 1E+7) {
            std::vector<double> c(boost::num_vertices(g));
            boost::brandes_betweenness_centrality(g, &c[0]);
            for (std::size_t i = 0; i != c.size(); ++i) {
                if (c[i] > maxCentrality) {
                    maxCentrality = c[i];
                    mostCentral = i;
                }
            }
        }
        bestRoot = mostCentral;
#else
        // Find the most extreme vertex.  This is a vertex that will almost
        // certainly be at the end of one or another track.
        std::size_t extremeVertex = 0;
        CP::THandle<CP::TReconCluster> extremeCluster;
        for (boost::tie(vi,vi_end) = boost::vertices(g); vi != vi_end; ++vi) {
            if (extremeCluster) {
                if (extremeCluster->GetPosition().X()
                    < g[*vi].cluster->GetPosition().X()) continue;
            }
            int edges = 0;
            MST::adjacency_iterator_t ai, ai_end;
            for (boost::tie(ai,ai_end) = boost::adjacent_vertices(*vi,g);
                 ai != ai_end; ++ai) {
                ++edges;
            }
            if (edges < 1) continue;
            extremeCluster = g[*vi].cluster;
            extremeVertex = *vi;
        }
        bestRoot = extremeVertex;
#endif

        // This finds the MST with the vertex bestRoot as the root, and uses
        // the EdgeProperties::length field for the edge length (the
        // weight_map needed by prim_minimum_spanning_tree).
        std::vector<MST::vertex_t> p(boost::num_vertices(g));
        boost::prim_minimum_spanning_tree(
            g, &p[0],
            boost::weight_map(boost::get(&MST::edge_properties_t::length,g)).
            distance_map(boost::get(&MST::vertex_properties_t::distance,g)).
            root_vertex(bestRoot));

        // Find the number of children, the number of steps to the root, the
        // distance to root, and the parent for each vertex.
        double maxDistance = 0.0;
        for (std::size_t i = 0; i != p.size(); ++i) {
            if (p[i] == i && g[i].distance>0) continue;
            g[i].parent = p[i];
            ++g[p[i]].children;
            g[i].rootDistance = g[i].distance;
            std::size_t j = i;
            while (j != p[j]) {
                g[i].rootDistance += g[p[j]].distance;
                ++g[i].steps;
                j = p[j];
            }
            maxDistance = std::max(maxDistance, g[i].rootDistance);
        }

        // Find the terminal vertices
        typedef std::pair<double,MST::vertex_t> termDist;
        std::vector<termDist> terminals;
        
        for (boost::tie(vi,vi_end) = boost::vertices(g); vi != vi_end; ++vi) {
            if (g[*vi].parent < 0) continue;
            if (g[*vi].children != 0) continue;
            terminals.push_back(std::make_pair(g[*vi].rootDistance, *vi));
        }
        std::sort(terminals.begin(), terminals.end());

        // Build tracks starting from the terminal vertex furthest from the
        // root.
        for (std::vector<termDist>::reverse_iterator t = terminals.rbegin();
             t != terminals.rend(); ++t) {
            std::vector< CP::THandle<CP::TReconBase> > nodes;
            int current = t->second; 
            for (;;) {
                // Check if we're at the root, or not connected.
                if (current < 0) break;
                CP::THandle<CP::TReconBase> object = g[current].cluster ;
                // Add the object to the vector of nodes.
                nodes.push_back(object);
                // The current vertex was already added to a track, so it's
                // the last one for the current track.
                if (g[current].marked) break;
                g[current].marked = true;
                // The current vertex is the root (and so is it's own parent).
                if (current == g[current].parent) break;
                current = g[current].parent;
            }

            // There are enough nodes in the new track to save.
            if (nodes.size() > 1) {
                CP::THandle<CP::TReconTrack> track 
                    = CreateTrackFromClusters("TMinimalSpanningTrack",
                                  nodes.begin(), nodes.end());
                final->push_back(track);
            }

        }

        // Move any remaining vertices out of the graph into
        // remainingClusters.
        remainingClusters.clear();
        for (boost::tie(vi,vi_end) = boost::vertices(g); vi != vi_end; ++vi) {
            if (g[*vi].marked) continue;
            remainingClusters.push_back(g[*vi].cluster);
        }

    } while ((0 < --throttle) && (2 < remainingClusters.size()));

    CaptNamedLog("MST","Clusters not added to any track: " 
                 << remainingClusters.size());

    // Save any remaining clusters into the output container.
    std::copy(remainingClusters.begin(), remainingClusters.end(), 
              std::back_inserter(*final));

    result->AddResultsContainer(final.release());

    return result;
}
