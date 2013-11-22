#include "TSplitTracks.hxx"
#include "TSegmentTrackFit.hxx"
#include "TTrackFit.hxx"
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

#include <memory>
#include <cmath>

CP::TSplitTracks::TSplitTracks()
    : TAlgorithm("TSplitTracks", 
                 "Split Tracks with Kinks and Gaps") {
    fThreeInLineCut = 3.0;
    fSplitDistanceCut = 50.0*unit::mm;
    fEndDistanceCut = 10.0*unit::mm;
}

CP::TSplitTracks::~TSplitTracks() { }

void CP::TSplitTracks::SaveTrack(
    CP::TReconObjectContainer& container,
    ClusterContainer::iterator begin, 
    ClusterContainer::iterator end) {

    // First make sure there are enough clusters for a track.
    if (end-begin < 2) {
        CaptNamedLog("Split", "Save track clusters.  To short");
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
            CaptNamedLog("Split", "Save track clusters."
                         "  Short track with gap of " << maxDist);
            std::copy(begin,end,std::back_inserter(container));
            return;
        }
    }

    CaptNamedLog("Split", "Save Track " << end-begin);

    CP::THandle<CP::TReconTrack> track
        = CP::CreateTrack("TSplitTracks",begin, end);
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

        CaptNamedLog("Split",
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
            CaptNamedLog("Split", "At front X2: " << v << "   D: " << d);
            if (v < fThreeInLineCut && d < fEndDistanceCut) break;
            final->push_back(*begin);
            ++begin;
            CaptNamedLog("Split",
                         "Discard a front cluster. New length: "<<end-begin);
        }

        // Check to see if the back clusters should be removed from the track.
        // Add any clusters that are removed to the final objects.
        while (end != begin+3) {
            double v = ThreeInLine(*(end-1), *(end-2), *(end-3));
            if (v < fThreeInLineCut) break;
            double d = ClusterDistance(**(end-1), **(end-2));
            CaptNamedLog("Split", "At back X2: " << v << "   D: " << d);
            if (v < fThreeInLineCut && d < fEndDistanceCut) break;
            final->push_back(*(end-1));
            --end;
            CaptNamedLog("Split",
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
                CaptNamedLog("Split","Discard a short track.  Length: "
                        << end-begin);
            }
            else if (end-begin == 3) {
                // These are three clusters in a row, so check if they are
                // consistent with a line.  If they are, make a track from the
                // clusters and save it to the final objects.
                double v = ThreeInLine(*begin, *(begin+1), *(begin+2));
                double d = std::max(ClusterDistance(**begin, **(begin+1)),
                                    ClusterDistance(**(begin+1), **(begin+2)));
                CaptNamedLog("Split", "Short check " << v << " " << d);
                if (v < fThreeInLineCut && d < fEndDistanceCut) {
                    SaveTrack(*final,begin,end);
                } 
                else {
                    std::copy(begin,end,std::back_inserter(*final));
                    CaptNamedLog("Split","Discard a short track.  Length: "
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
                CaptNamedLog("Split", "Medium check " << v << " " << d);
                if (v < fThreeInLineCut && d < fEndDistanceCut) {
                    SaveTrack(*final,begin,end);
                } 
                else {
                    std::copy(begin,end,std::back_inserter(*final));
                    CaptNamedLog("Split","Discard a short track.  Length: "
                            << end-begin);
                }
            }
            continue;
        }

        // Find the sharpest kink in the track.  There are at least five
        // clusters in the track when this bit of code starts.  The kink can't
        // be right at either end.
        ClusterContainer::iterator sharpestKink;
        double worstChi2 = 0.0;
        for (ClusterContainer::iterator i = begin+1; i+3 != end; ++i) {
            double v = ThreeInLine(*(i),*(i+1), *(i+2));
            if (worstChi2 < v) {
                worstChi2 = v;
                sharpestKink = i+1;
            }
        }
        CaptNamedLog("Split", "Three in line for split: " << worstChi2);

        // If the worstChi2 is too big, there is a kink and the track should
        // be split.  The two pieces are put back on the stack and we start
        // over again.
        if (worstChi2 > fThreeInLineCut) {
            track = CreateTrack("TSplitTracks",begin,sharpestKink+1);
            trackStack.push_back(track);
            track = CreateTrack("TSplitTracks",sharpestKink,end);
            trackStack.push_back(track);
            CaptNamedLog("Split","Split a track."
                         << " X2: " << worstChi2
                         << "  Original: "   << end-begin
                         << "  First: "   << sharpestKink-begin+1
                         << "  Second: "   << end-sharpestKink);
            continue;
        }

        // Find the biggest gap in the track.  There are at least five
        // clusters in the track when this bit of code starts.  
        ClusterContainer::iterator biggestGap;
        double maxDist = 0.0;
        for (ClusterContainer::iterator i = begin+1; i+3 != end; ++i) {
            double v = ClusterDistance(**(i),**(i+1));
            if (maxDist < v) {
                maxDist = v;
                biggestGap = i+1;
            }
        }
        CaptNamedLog("Split", "Cluster distance for split: " << maxDist);

        // If the maxDist is too big, there is a kink and the track should
        // be split.  The two pieces are put back on the stack and we start
        // over again.
        if (maxDist > fSplitDistanceCut) {
            track = CreateTrack("TSplitTracks",begin,biggestGap);
            trackStack.push_back(track);
            track = CreateTrack("TSplitTracks",biggestGap,end);
            trackStack.push_back(track);
            CaptNamedLog("Split","Split a track at a gap."
                         << "  Gap: " << maxDist 
                         << "  Original: "   << end-begin
                         << "  First: "   << biggestGap-begin
                         << "  Second: "   << end-biggestGap);
            continue;
        }
        
        // Whatever gets here should just be saved...
        SaveTrack(*final,begin,end);
    }
    
    result->AddResultsContainer(final.release());

    return result;
}
