#include "TBestTube.hxx"
#include "TTubePredicate.hxx"
#include "TMajorAxisComparator.hxx"
#include "ostreamTVector3.hxx"

#include "TCaptLog.hxx"

#include <TRandom.h>

CP::TReconObjectContainer::iterator 
CP::TBestTube::Randomize(CP::TReconObjectContainer::iterator i, 
                         CP::TReconObjectContainer::iterator begin,
                         CP::TReconObjectContainer::iterator end) {
    int diff = std::distance(begin,end);
    int step = gRandom->Integer(diff);
    i = begin;
    std::advance(i,step);
    return i;
}

double CP::TBestTube::TubeWeight(const TVector3& end1, const TVector3& end2,
                                 CP::TReconObjectContainer::iterator begin, 
                                 CP::TReconObjectContainer::iterator end) {
    TTubePredicate tube(end1,end2);
    double result = 0.0;
    while (begin != end) {
        CP::THandle<CP::TReconCluster> b(*begin);
        if (tube(b->GetPosition().Vect())) {
            result += 1.0;
        }
        ++begin;
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
    CP::TReconObjectContainer::iterator best1;
    CP::TReconObjectContainer::iterator best2;
    double bestWeight = 0.0;
    
    int diff = std::distance(begin,end);
    CP::TReconObjectContainer::iterator end1 = begin;
    std::advance(end1,diff/3);
    CP::TReconObjectContainer::iterator end2 = begin; 
    std::advance(end2,2*diff/3);
    for (int trial=0; trial<diff; ++trial) {
        double weight = TubeWeight(end1,end2,begin,end);
        if (weight > bestWeight) {
            best1 = end1;
            best2 = end2;
            bestWeight = weight;
        }
        end1 = Randomize(end1,begin,end);
        end2 = Randomize(end2,begin,end);
    }

    CP::THandle<CP::TReconCluster> tmp = *end1;
    fEnd1 = tmp->GetPosition().Vect();
    tmp = *end2;
    fEnd2 = tmp->GetPosition().Vect();

    CaptNamedDebug("tube","Found seed with " << bestWeight 
               << " from " << fEnd1 << "-->" << fEnd2);

}

void CP::TBestTube::FillSeed(CP::TReconObjectContainer& seed) {
    seed.clear();
    TTubePredicate tube(fEnd1,fEnd2);
    for (CP::TReconObjectContainer::iterator c = fClusterContainer.begin();
         c != fClusterContainer.end(); ++c) {
        CP::THandle<CP::TReconCluster> cluster(*c);
        if (tube(cluster->GetPosition().Vect())) {
            seed.push_back(cluster);
        }
    }
    CaptNamedDebug("tube","Seed with " << seed.size() << " clusters");

    CP::TMajorAxisComparator axis(seed);
    std::sort(seed.begin(), seed.end(), axis);

    CP::THandle<CP::TReconCluster> c1(seed.front());
    CP::THandle<CP::TReconCluster> c2(seed.back());
    CaptNamedDebug("tube","   from " << c1->GetPosition().Vect());
    CaptNamedDebug("tube","     to " << c2->GetPosition().Vect());
}

void CP::TBestTube::FillRemains(CP::TReconObjectContainer& remains) {
    remains.clear();
    TTubePredicate tube(fEnd1,fEnd2);
    for (CP::TReconObjectContainer::iterator c = fClusterContainer.begin();
         c != fClusterContainer.end(); ++c) {
        CP::THandle<CP::TReconCluster> cluster(*c);
        if (! tube(cluster->GetPosition().Vect())) {
            remains.push_back(cluster);
        }
    }
    CaptNamedDebug("tube","Remaining clusters " << remains.size());
}

