#include "TBestTube.hxx"
#include "TTubePredicate.hxx"
#include "TMajorAxisComparator.hxx"
#include "ostreamTVector3.hxx"

#include "TCaptLog.hxx"

namespace {
    // This finds the minimum distance between hits in the two clusters.
    double ClusterDistance(const CP::TReconCluster& a, 
                           const CP::TReconCluster& b) {
        CP::THandle<CP::THitSelection> aHits = a.GetHits();
        CP::THandle<CP::THitSelection> bHits = b.GetHits();
        double minDist = 1000*unit::meter;
        for (CP::THitSelection::iterator j = aHits->begin(); 
             j != aHits->end(); ++j) {
            for (CP::THitSelection::iterator k = bHits->begin(); 
                 k != bHits->end(); ++k) {
                if (j == k) continue;
                double dist = ((*j)->GetPosition()-(*k)->GetPosition()).Mag();
                if (dist < minDist) minDist = dist;
            }
        }
        return minDist;
    }

    // Sort the clusters in order from end1 to end2.
    struct TubeSort {
    public:
        TubeSort(TVector3 end1, TVector3 end2) {
            fStart = end1;
            fDir = (end2-end1).Unit();
        }

        bool operator () (const CP::THandle<CP::TReconCluster>& lhs,
                          const CP::THandle<CP::TReconCluster>& rhs) {
            TVector3 lpos = lhs->GetPosition().Vect();
            double lval = (lpos-fStart)*fDir;

            TVector3 rpos = rhs->GetPosition().Vect();
            double rval = (rpos-fStart)*fDir;
            return lval < rval;
        }

    private:
        TVector3 fStart;
        TVector3 fDir;
    };
}

void CP::TBestTube::FindTube(const TVector3& end1, const TVector3& end2,
                             CP::TReconObjectContainer::iterator begin, 
                             CP::TReconObjectContainer::iterator end,
                             CP::TReconObjectContainer& output) {
    TTubePredicate tube(end1,end2);

    std::vector< CP::THandle<CP::TReconCluster> > clusters;
    while (begin != end) {
        CP::THandle<CP::TReconCluster> b(*begin);
        if (tube(b->GetPosition().Vect())) {
            output.push_back(b);
        }
        ++begin;
    }

    std::sort(output.begin(), output.end(), TubeSort(end1,end2));
}

double CP::TBestTube::TubeWeight(const TVector3& end1, const TVector3& end2,
                                 CP::TReconObjectContainer::iterator begin, 
                                 CP::TReconObjectContainer::iterator end) {
    TReconObjectContainer tube;
    tube.reserve(50);

    // Find the clusters that are in the tube.
    FindTube(end1, end2, begin, end, tube);

    // Figure out if the hits in the tube make a good seed.
    if (tube.size() < 3) return 0.0;

    // Check that the hits are all neighbors.
    for (CP::TReconObjectContainer::iterator c = tube.begin()+1;
         c != tube.end(); ++c) {
        CP::TReconObjectContainer::iterator b = c-1;
        CP::THandle<CP::TReconCluster> f(*b);
        CP::THandle<CP::TReconCluster> s(*c);
        // Oops.  This isn't a cluster.  This can't be a seed.
        if (!f || !s) {
            CaptError("TBestTube called with non-cluster in container.");
            return 0.0;
        }
        double dist = ClusterDistance(*f, *s);
        // Make sure the seed is continuous (i.e. no "big" gaps).
        if (dist > 5*unit::mm) return 0.0;
    }

    // We have a good seed candidate, so find it's weight.
    double result = tube.size();

    // We don't want the seed to be too long...
    if (result > 10.0) {
        result = 10.0*exp(-0.1*(result-10.0));
    }

    return result;
}

double CP::TBestTube::TubeWeight(CP::TReconObjectContainer::iterator end1,
                                 CP::TReconObjectContainer::iterator end2,
                                 CP::TReconObjectContainer::iterator begin, 
                                 CP::TReconObjectContainer::iterator end) {
    if (end1 == end2) return 0.0;
    CP::THandle<CP::TReconCluster> e1(*end1);
    CP::THandle<CP::TReconCluster> e2(*end2);
    return TubeWeight(e1->GetPosition().Vect(), e2->GetPosition().Vect(),
                      begin, end);
}

void CP::TBestTube::Process(const CP::TReconObjectContainer& input) {
    fClusterContainer.clear();
    fClusterContainer.reserve(input.size());
    std::copy(input.begin(),input.end(), std::back_inserter(fClusterContainer));
    CP::TReconObjectContainer::iterator begin = fClusterContainer.begin();
    CP::TReconObjectContainer::iterator end = fClusterContainer.end();
    CP::TReconObjectContainer::iterator best1 = end;
    CP::TReconObjectContainer::iterator best2 = end;
    double bestWeight = 0.0;
    fFoundTube = false;

    int diff = std::distance(begin,end);
    CaptNamedDebug("bestTube","Input clusters " << diff);
    CP::TReconObjectContainer::iterator end1 = begin;
    CP::TReconObjectContainer::iterator end2 = begin;
    std::advance(end2,2*diff/3);
    double trials = std::pow(1.0*diff,1.5) + 1;
    CaptNamedDebug("bestTube","Trials " << int(trials));
    for (CP::TReconObjectContainer::iterator end1 = begin;
         end1 != end;  ++end1) {
        for (CP::TReconObjectContainer::iterator end2 = end1+1;
             end2 != end; ++end2) {
            double weight = TubeWeight(end1,end2,begin,end);
            if (weight > bestWeight) {
                CP::THandle<CP::TReconCluster> t1 = *end1;
                CP::THandle<CP::TReconCluster> t2 = *end2;
                TVector3 p1 = t1->GetPosition().Vect();
                TVector3 p2 = t2->GetPosition().Vect();
                CaptNamedDebug("bestTube","New Weight "
                             << weight
                             << " from " << p1
                             << "-->" << p2
                             << " length is " << (p1-p2).Mag());
                best1 = end1;
                best2 = end2;
                fFoundTube = true;
                bestWeight = weight;
            }
        }
    }

    if (!fFoundTube) {
        CaptNamedInfo("bestTube","No seed found");
        return;
    }

    CP::THandle<CP::TReconCluster> tmp = *best1;
    fEnd1 = tmp->GetPosition().Vect();
    tmp = *best2;
    fEnd2 = tmp->GetPosition().Vect();

    CaptNamedInfo("bestTube","Found seed with " << bestWeight 
               << " from " << fEnd1 << "-->" << fEnd2);

}

void CP::TBestTube::FillSeed(CP::TReconObjectContainer& seed) {
    seed.clear();
    if (!fFoundTube) return;
    CaptNamedDebug("bestTube","Find hits for " << fEnd1 << "-->" << fEnd2);
    TTubePredicate tube(fEnd1,fEnd2);
    for (CP::TReconObjectContainer::iterator c = fClusterContainer.begin();
         c != fClusterContainer.end(); ++c) {
        CP::THandle<CP::TReconCluster> cluster(*c);
        if (tube(cluster->GetPosition().Vect())) {
            seed.push_back(cluster);
        }
    }
    CaptNamedDebug("bestTube","Seed with " << seed.size() << " clusters");

    CP::TMajorAxisComparator axis(seed);
    std::sort(seed.begin(), seed.end(), axis);

    CP::THandle<CP::TReconCluster> c1(seed.front());
    CP::THandle<CP::TReconCluster> c2(seed.back());
}

void CP::TBestTube::FillRemains(CP::TReconObjectContainer& remains) {
    remains.clear();
    if (!fFoundTube) {
        std::copy(fClusterContainer.begin(), fClusterContainer.end(),
                  std::back_inserter(remains));
        return;
    }
    TTubePredicate tube(fEnd1,fEnd2);
    for (CP::TReconObjectContainer::iterator c = fClusterContainer.begin();
         c != fClusterContainer.end(); ++c) {
        CP::THandle<CP::TReconCluster> cluster(*c);
        if (! tube(cluster->GetPosition().Vect())) {
            remains.push_back(cluster);
        }
    }
    CaptNamedDebug("bestTube","Remaining clusters " << remains.size());
}

