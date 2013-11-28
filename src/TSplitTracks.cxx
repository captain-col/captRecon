#include "TSplitTracks.hxx"
#include "TTrackFit.hxx"
#include "TSegmentTrackFit.hxx"
#include "ClusterDistance.hxx"
#include "CreateTrack.hxx"

#include <THandle.hxx>
#include <TReconTrack.hxx>
#include <TReconCluster.hxx>
#include <TCaptLog.hxx>
#include <HEPUnits.hxx>
#include <TUnitsTable.hxx>
#include <TRuntimeParameters.hxx>

#include <TMatrixD.h>
#include <TPrincipal.h>

#include <memory>
#include <cmath>

namespace {
    struct CompareRecon {
        bool operator () (CP::THandle<CP::TReconBase> lhs, 
                          CP::THandle<CP::TReconBase> rhs) {
            CP::THandle<CP::TReconTrack> lt = lhs;
            CP::THandle<CP::TReconTrack> rt = rhs;
            CP::THandle<CP::TReconCluster> lc = lhs;
            CP::THandle<CP::TReconCluster> rc = rhs;
            
            if (lt && rt) {
                return lt->GetNodes().size() > rt->GetNodes().size();
            }
            
            if (lt && !rt) return true;
            if (!lt && rt) return false;

            if (lc && rc) {
                return lc->GetEDeposit() > rc->GetEDeposit();
            }                
            
            if (lc && !rc) return true;
            if (!lc && rc) return false;

            return CP::GetPointer(lhs) < CP::GetPointer(rhs);
        }
    };
};

CP::TSplitTracks::TSplitTracks()
    : TAlgorithm("TSplitTracks", 
                 "Split Tracks with Kinks and Gaps") {
    fThreeInLineCut = 3.0;
    fSplitDistanceCut = 30.0*unit::mm;
    fEndDistanceCut = 10.0*unit::mm;
    fKinkAngleCut = 20*unit::degree;
}

CP::TSplitTracks::~TSplitTracks() { }

void CP::TSplitTracks::SaveTrack(
    CP::TReconObjectContainer& container,
    ClusterContainer::iterator begin, 
    ClusterContainer::iterator end) {

    // First make sure there are enough clusters for a track.
    if (end-begin < 2) {
        CaptNamedInfo("Split", "Save track clusters.  To short");
        std::copy(begin,end,std::back_inserter(container));
        return;
    }

    // There are enough clusters, but for a short track make sure they are "in
    // contact".
    if (end-begin < 5) {
        double maxDist = 0.0;
        for (ClusterContainer::iterator i = begin; i+1 != end; ++i) {
            maxDist = std::max(maxDist, CP::ClusterDistance(**(i), **(i+1)));
        }
        if (maxDist > fSplitDistanceCut) {
            CaptNamedInfo("Split", "Save track clusters."
                         "  Short track with gap of " << maxDist);
            std::copy(begin,end,std::back_inserter(container));
            return;
        }
    }

    CP::THandle<CP::TReconTrack> track
        = CP::CreateTrack("TSplitTracks",begin, end);

    if (!track) {
        CaptNamedLog("Split", "Track not created");
    }

    CP::TTrackFit fitter;
    track = fitter(track);

    container.push_back(track);
}

double
CP::TSplitTracks::ThreeInLine(CP::THandle<CP::TReconCluster> a,
                              CP::THandle<CP::TReconCluster> b, 
                              CP::THandle<CP::TReconCluster> c) {
    // The assumption is that B lies between A and C, but the code will work
    // even if that's not true.

    // Find the deviation between the b cluster and the line between the a and
    // c clusters.
    TVector3 diffBA = b->GetPosition().Vect()-a->GetPosition().Vect();
    TVector3 dirCA = (c->GetPosition().Vect()-a->GetPosition().Vect()).Unit();
    double proj = diffBA*dirCA;
    TVector3 deltaB = diffBA - proj*dirCA;

    // Find the covariance.
    TMatrixD cov(3,3);
    for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
            cov(i,j) = 0.0;
            cov(i,j) += a->GetState()->GetPositionCovariance(i,j);
            cov(i,j) += b->GetState()->GetPositionCovariance(i,j);
            cov(i,j) += c->GetState()->GetPositionCovariance(i,j);
        }
    }

    // Invert to find the error matrix.  This modifies the covariance matrix.
    cov.InvertFast();
    
    // Find the chi2...
    double result = 0.0;
    for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
            result += deltaB[i]*cov(i,j)*deltaB[j];
        }
    }

    return result;
}

double CP::TSplitTracks::KinkAngle(ClusterContainer::iterator here, 
                                   ClusterContainer::iterator begin,
                                   ClusterContainer::iterator end) {
    // The minimum and maximum number of clusters to include in the segments
    // to either side of the cluster being checked.  The count includes the
    // current cluster.
    int minStep = 5;
    int maxStep = 10;

    // No kinks right at the end.
    if (here-begin+1 < minStep) return 0;
    if (end-here < minStep) return 0;

    double p0[3] = {0.0,0,0};
    double p1[3] = {1.0,0,0};
    double x[3];

    // Do a simple fit in the backward direction.  This includes the cluster
    // being checked.
    int backStep = here-begin;
    if (backStep > maxStep) backStep = maxStep;
    TPrincipal pca1(3,"");
    for (ClusterContainer::iterator i = here-backStep; i != here+1; ++i) {
        double row[3] 
            = {(*i)->GetPosition().X(),
               (*i)->GetPosition().Y(),
               (*i)->GetPosition().Z()};
        pca1.AddRow(row);
    }
    pca1.MakePrincipals();
    pca1.P2X(p0,x,1);
    TVector3 base1(x);
    pca1.P2X(p1,x,1);
    TVector3 dir1(x); 
    dir1 = (dir1-base1).Unit();

    // Do a simple fit in the forward direction.  This includes the cluster
    // being checked.
    int foreStep = end-here;
    if (foreStep > maxStep) foreStep = maxStep;
    TPrincipal pca2(3,"");
    for (ClusterContainer::iterator i = here; i != here+foreStep; ++i) {
        double row[3] 
            = {(*i)->GetPosition().X(),
               (*i)->GetPosition().Y(),
               (*i)->GetPosition().Z()};
        pca2.AddRow(row);
    }
    pca2.MakePrincipals();
    pca2.P2X(p0,x,1);
    TVector3 base2(x);
    pca2.P2X(p1,x,1);
    TVector3 dir2(x); 
    dir2 = (dir2-base2).Unit();

    // Find the cosine of the angle between the fore and back segments.
    double dcos = std::abs(dir1*dir2);
    if (dcos >= 1.0) return 0.0;

    // Return the angle between the segments.  This could be done with the
    // sine (using a sqrt), but I think acos is fast enough.
    return std::acos(dcos);
}

CP::THandle<CP::TAlgorithmResult>
CP::TSplitTracks::Process(const CP::TAlgorithmResult& input,
                          const CP::TAlgorithmResult&,
                          const CP::TAlgorithmResult&) {
    
    CP::THandle<CP::TReconObjectContainer> inputObjects 
        = input.GetResultsContainer();

    CaptLog("TSplitTracks Process " << GetEvent().GetContext());

    if (!inputObjects) {
        CaptError("No input objects");
        return CP::THandle<CP::TAlgorithmResult>();
    }

    // Create the output containers.
    CP::THandle<CP::TAlgorithmResult> result = CreateResult();
    std::auto_ptr<CP::TReconObjectContainer> 
        final(new CP::TReconObjectContainer("final"));

    // Create a stack to keep tracks in.  When the stack is empty, all the
    // tracks that need to be split have been split.
    typedef std::vector< CP::THandle<CP::TReconTrack> > TrackStack;
    TrackStack trackStack(inputObjects->size());
    trackStack.clear();

    // Make a copy of all of the tracks in the input objects, and save the
    // non-tracks to the output.
    for (CP::TReconObjectContainer::iterator t = inputObjects->begin();
         t != inputObjects->end(); ++t) {
        CP::THandle<CP::TReconTrack> track = *t;
        if (!track) {
            final->push_back(*t);
            continue;
        }
        trackStack.push_back(track);
    }

    // Pop a track off the stack and see if it should be split.  If the track
    // is split, the two halfs are both pushed back on the stack.  If it
    // doesn't get split, the track get's pushed into the final object
    // container.
    while (!trackStack.empty()) {
        CP::THandle<CP::TReconTrack> track = trackStack.back();
        trackStack.pop_back();

        CaptNamedInfo("Split",
                     "Track Stack: " << trackStack.size()
                     << "    Track size: " << track->GetNodes().size());

        // Make a copy of the clusters in the track.  These are basically the
        // track nodes.
        ClusterContainer clusters;
        for (CP::TReconNodeContainer::iterator n = track->GetNodes().begin();
             n != track->GetNodes().end(); ++n) {
            CP::THandle<CP::TReconCluster> c = (*n)->GetObject();
            if (!c) {
                CaptError("Track node not made from a cluster");
                throw;
            }
            clusters.push_back(c);
        }
            
        ClusterContainer::iterator begin, end;
        begin = clusters.begin();
        end = clusters.end();

        // Protect against small tracks...
        if (end-begin < 3) {
            SaveTrack(*final,begin,end);
            continue;
        }

        // Check to see if the front clusters should be removed from the
        // track.  Add any clusters that are removed to the final objects.
        while ((begin + 3) != end) {
            double v = ThreeInLine(*begin, *(begin+1), *(begin+2));
            double d = ClusterDistance(**begin, **(begin+1));
            CaptNamedInfo("Split", "At front X2: " << v << "   D: " << d);
            if (v < fThreeInLineCut && d < fEndDistanceCut) break;
            final->push_back(*begin);
            ++begin;
            CaptNamedInfo("Split",
                         "Discard a front cluster. New length: "<<end-begin);
        }

        // Check to see if the back clusters should be removed from the track.
        // Add any clusters that are removed to the final objects.
        while (end != begin+3) {
            double v = ThreeInLine(*(end-1), *(end-2), *(end-3));
            if (v < fThreeInLineCut) break;
            double d = ClusterDistance(**(end-1), **(end-2));
            CaptNamedInfo("Split", "At back X2: " << v << "   D: " << d);
            if (v < fThreeInLineCut && d < fEndDistanceCut) break;
            final->push_back(*(end-1));
            --end;
            CaptNamedInfo("Split",
                         "Discard a back cluster. New length: "<<end-begin);
        }

        // Check to see if the track length has been reduced so much that it
        // should be split up into it's clusters.  The track will either be
        // added to the final set of objects, or split into clusters which are
        // added to the final set of objects.
        if (end-begin < 5) {
            if (end-begin < 3) {
                // This is a very short track stub after removing clusters, so
                // just save the clusters.
                std::copy(begin,end,std::back_inserter(*final));
                CaptNamedInfo("Split","Discard a short track.  Length: "
                        << end-begin);
            }
            else if (end-begin == 3) {
                // These are three clusters in a row, so check if they are
                // consistent with a line.  If they are, make a track from the
                // clusters and save it to the final objects.
                double v = ThreeInLine(*begin, *(begin+1), *(begin+2));
                double d = std::max(ClusterDistance(**begin, **(begin+1)),
                                    ClusterDistance(**(begin+1), **(begin+2)));
                CaptNamedInfo("Split", "Short check " << v << " " << d);
                if (v < fThreeInLineCut && d < fEndDistanceCut) {
                    SaveTrack(*final,begin,end);
                } 
                else {
                    std::copy(begin,end,std::back_inserter(*final));
                    CaptNamedInfo("Split","Discard a short track.  Length: "
                            << end-begin);
                }
            }
            else if (end-begin == 4) {
                // These are four clusters in a row, so check if they are
                // consistent with a line.  If they are, make a track from the
                // clusters and save it to the final objects.
                double v1 = ThreeInLine(*begin, *(begin+1), *(begin+2));
                double v2 = ThreeInLine(*(begin+1), *(begin+2), *(begin+3));
                double v = std::max(v1,v2);
                double d = std::max(ClusterDistance(**begin, **(begin+1)),
                                    ClusterDistance(**(begin+2), **(begin+3)));
                CaptNamedInfo("Split", "Medium check " << v << " " << d);
                if (v < fThreeInLineCut && d < fEndDistanceCut) {
                    SaveTrack(*final,begin,end);
                } 
                else {
                    std::copy(begin,end,std::back_inserter(*final));
                    CaptNamedInfo("Split","Discard a short track.  Length: "
                            << end-begin);
                }
            }
            continue;
        }

        // Look for parts of the track where three clusters are not in a line.
        // The segment with the worst chi2 will be a "short" kink in the
        // track.  There are at least five clusters in the track when this bit
        // of code starts.  The kink can't be right at either end.
        ClusterContainer::iterator sharpestKink;
        double worstChi2 = 0.0;
        for (ClusterContainer::iterator i = begin+1; i+3 != end; ++i) {
            double v = ThreeInLine(*(i),*(i+1), *(i+2));
            if (worstChi2 < v) {
                worstChi2 = v;
                sharpestKink = i+1;
            }
        }
        CaptNamedInfo("Split", "Three in line for split: " << worstChi2);

        // If the worstChi2 is too big, there is a kink and the track should
        // be split.  The two pieces are put back on the stack and we start
        // over again.
        if (worstChi2 > fThreeInLineCut) {
            track = CreateTrack("TSplitTracks",begin,sharpestKink+1);
            trackStack.push_back(track);
            track = CreateTrack("TSplitTracks",sharpestKink,end);
            trackStack.push_back(track);
            CaptNamedInfo("Split","Split a track."
                         << " X2: " << worstChi2
                         << "  Original: "   << end-begin
                         << "  First: "   << sharpestKink-begin+1
                         << "  Second: "   << end-sharpestKink);
            continue;
        }

        // Find the biggest gap in the track.  The can't be right at the end
        // of the track.  There are at least 5 clusters in the track to get to
        // this point, so the check (i+3 != end) should be OK.
        ClusterContainer::iterator biggestGap;
        double maxDist = 0.0;
        for (ClusterContainer::iterator i = begin+1; i+3 != end; ++i) {
            double v = ClusterDistance(**(i),**(i+1));
            if (maxDist < v) {
                maxDist = v;
                biggestGap = i+1;
            }
        }
        CaptNamedInfo("Split", "Cluster distance for split: " << maxDist);

        // If the biggest gap is too big the track should be split.  The two
        // pieces are put back on the stack and we start over again.
        if (maxDist > fSplitDistanceCut) {
            track = CreateTrack("TSplitTracks",begin,biggestGap);
            trackStack.push_back(track);
            track = CreateTrack("TSplitTracks",biggestGap,end);
            trackStack.push_back(track);
            CaptNamedInfo("Split","Split a track at a gap."
                         << "  Gap: " << maxDist 
                         << "  Original: "   << end-begin
                         << "  First: "   << biggestGap-begin
                         << "  Second: "   << end-biggestGap);
            continue;
        }
        CaptNamedInfo("Split", "Cluster distance for split: " << maxDist);

        // Look through the clusters and find the one with the sharpest kink
        // based on the track segments to either side of the current cluster.
        // This reuses the sharpestKink iterator.
        double biggestAngle = 0.0;
        for (ClusterContainer::iterator i = begin; i != end; ++i) {
            double v = KinkAngle(i,begin,end);
            if (biggestAngle < v) {
                biggestAngle = v;
                sharpestKink = i;
            }
        }
        CaptNamedLog("Split","Split a track."
                     << " Angle: " << biggestAngle
                     << "  max angle: " << fKinkAngleCut);

        // If the kink angle is to big there is a kink and the track should
        // be split.  The two pieces are put back on the stack and we start
        // over again.
        if (biggestAngle > fKinkAngleCut) {
            track = CreateTrack("TSplitTracks",begin,sharpestKink+1);
            trackStack.push_back(track);
            track = CreateTrack("TSplitTracks",sharpestKink,end);
            trackStack.push_back(track);
            CaptNamedLog("Split","Split a track."
                         << " Angle: " << biggestAngle
                         << "  Original: "   << end-begin
                         << "  First: "   << sharpestKink-begin+1
                         << "  Second: "   << end-sharpestKink);
            continue;
        }

        // Whatever gets here should just be saved...
        SaveTrack(*final,begin,end);
    }
    
    CaptError("Sort the tracks");

    std::sort(final->begin(), final->end(), CompareRecon());

    result->AddResultsContainer(final.release());

    return result;
}
