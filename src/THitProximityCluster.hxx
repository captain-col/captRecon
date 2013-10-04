#ifndef THitProximityCluster_hxx_seen
#include "TTmplDensityCluster.hxx"

#include <THit.hxx>
#include <THandle.hxx>

namespace CP {
    namespace HitProximity {
        /// The type of argument accepted by the clustering.
        typedef CP::THandle<CP::THit> Arg;

        /// Check to see if two hits are close together.  Since the hits might
        /// have very different extents in the XY plane and along the Z axis,
        /// this subtracts off the hit sizes.
        struct Metric {
            double operator() (const Arg& lhs, const Arg& rhs);
        };

        /// The class to apply a density clustering to hits.
        typedef CP::TTmplDensityCluster< Arg, Metric > Cluster;
    };
};

#endif

