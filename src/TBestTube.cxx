#include "TBestTube.hxx"
#include "TTubePredicate.hxx"
#include "TMajorAxisComparator.hxx"
#include "ClusterDistance.hxx"

#include <TCaptLog.hxx>
#include <ostreamTVector3.hxx>

#include <TRandom.h>


namespace {
    // Choose a random iterator in the list of clusters.
    CP::TReconObjectContainer::iterator 
    Randomize(CP::TReconObjectContainer::iterator begin,
              CP::TReconObjectContainer::iterator end) {
        int diff = std::distance(begin,end);
        int step = gRandom->Integer(diff);
        return begin + step;
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

CP::TBestTube::TBestTube() 
    :  fGoodEnoughTube(5.0), fFoundTube(false), fMaxTrials(100*100) { }

CP::TBestTube::~TBestTube() { }

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
        double dist = CP::ClusterDistance(*f, *s);
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
    CaptNamedDebug("bestTube","Input clusters: " << diff);

    // Find the number of clusters to be checked, and calculate a maximum
    // number of trials.  For "small" events, this does an exhaustive search.
    // If the event gets to bit, this does a random search.
    int trials = diff*diff;
    bool limitTrials = false;
    if (trials > fMaxTrials) {
        // Limit trials to the number of clusters^(3/2).  This is the
        // recommended power in the literature for a random tube search in 3D.
        trials = std::pow(1.0*diff,1.5);
        // Make sure we do at least fMaxTrials checks.
        trials = std::max(trials,fMaxTrials);
        limitTrials = true;
        CaptNamedDebug("bestTube","Trials limited to " << trials);
    }


    int trial = 0; 
    bool terminateLoop = false;  // A poor man's try-catch...
    for (CP::TReconObjectContainer::iterator e1 = begin;
         e1 != end;  ++e1) {
        for (CP::TReconObjectContainer::iterator e2 = e1+1;
             e2 != end; ++e2) { 
            ++trial;
            // Make a copy of the loop iterators so they can be overridden.
            CP::TReconObjectContainer::iterator end1 = e1;
            CP::TReconObjectContainer::iterator end2 = e2;
            if (limitTrials) {
                end1 = Randomize(begin,end);
                end2 = Randomize(begin,end);
            }
            double weight = TubeWeight(end1,end2,begin,end);
            if (weight > bestWeight) {
                CP::THandle<CP::TReconCluster> t1 = *end1;
                CP::THandle<CP::TReconCluster> t2 = *end2;
                CaptNamedDebug("bestTube","New Weight "
                             << weight
                             << " from " << t1->GetPosition().Vect()
                             << "-->" << t2->GetPosition().Vect()
                             << " length is " 
                               << (t1->GetPosition().Vect()
                                   -t2->GetPosition().Vect()).Mag());
                best1 = end1;
                best2 = end2;
                fFoundTube = true;
                bestWeight = weight;
            }
            if (bestWeight > fGoodEnoughTube) terminateLoop = true;
            if (trials <= trial) terminateLoop = true;
            if (terminateLoop) break;
        }
        if (terminateLoop) break;
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

