#include "TMergeTracks.hxx"
#include "TTrackFit.hxx"
#include "ClusterDistance.hxx"
#include "CreateTrack.hxx"
#include "CompareReconObjects.hxx"
#include "TTrackMassFit.hxx"

#include <THandle.hxx>
#include <TReconTrack.hxx>
#include <TReconCluster.hxx>
#include <TCaptLog.hxx>
#include <HEPUnits.hxx>
#include <TUnitsTable.hxx>
#include <TRuntimeParameters.hxx>
#include <ostreamTVector3.hxx>
 
#include <TMatrixD.h>

#include <memory>
#include <deque>
#include <list>
#include <cmath>

CP::TMergeTracks::TMergeTracks()
    : TAlgorithm("TMergeTracks", 
                 "Merge Tracks with Kinks and Gaps") {
    fGoodnessCut = 16.0;
    fMergeDistanceCut = 5.0*unit::cm;
}

CP::TMergeTracks::~TMergeTracks() { }

CP::TMergeTracks::Orientation
CP::TMergeTracks::TrackOrientation(CP::THandle<CP::TReconTrack> first,
                                   CP::THandle<CP::TReconTrack> second) {
    /// The positions of the front and back of the tracks.
    TVector3 firstFront(first->GetFront()->GetPosition().Vect());
    TVector3 firstBack(first->GetBack()->GetPosition().Vect());
    TVector3 secondFront(second->GetFront()->GetPosition().Vect());
    TVector3 secondBack(second->GetBack()->GetPosition().Vect());
    
    Orientation orient 
        = TrackOrientation(firstFront,firstBack,secondFront,secondBack);

    if (orient != kNotClose) return orient;

    // The track ends were not close together, but the track ends are always
    // extrapolated to the edge of the cluster, so now check if the clusters
    // are close.  The objects *must* be clusters.
    CP::THandle<CP::TReconCluster> tempCluster;
    tempCluster = first->GetNodes().front()->GetObject();
    firstFront = tempCluster->GetPosition().Vect();

    tempCluster = first->GetNodes().back()->GetObject();
    firstBack = tempCluster->GetPosition().Vect();

    tempCluster = second->GetNodes().front()->GetObject();
    secondFront = tempCluster->GetPosition().Vect();

    tempCluster = second->GetNodes().back()->GetObject();
    secondBack = tempCluster->GetPosition().Vect();

    return TrackOrientation(firstFront,firstBack,secondFront,secondBack);
}

CP::TMergeTracks::Orientation
CP::TMergeTracks::TrackOrientation(const TVector3& firstFront,
                                   const TVector3& firstBack,
                                   const TVector3& secondFront,
                                   const TVector3& secondBack) {

    // Now figure out distances between the track ends.  The two closest
    // distance determine how the tracks will be arranged.
    double ffDist = (firstFront - secondFront).Mag();
    double fbDist = (firstFront - secondBack).Mag();
    double bfDist = (firstBack - secondFront).Mag();
    double bbDist = (firstBack - secondBack).Mag();

    double minDist = std::min(ffDist,std::min(fbDist,std::min(bfDist, bbDist)));
    if (minDist > fMergeDistanceCut) return kNotClose;

    // Figure out which ends are closest and use that to determine which order
    // the clusters get added.
    if (ffDist<fbDist && ffDist<bfDist && ffDist<bbDist) {
        // A front to front track.
        CaptNamedVerbose("Merge","Orient FF " << ffDist
                      << " " << fbDist
                      << " " << bfDist
                      << " " << bbDist);
        return kFrontFront;
    }
    else if (fbDist<ffDist && fbDist<bfDist && fbDist<bbDist) {
        // A front to back track.
        CaptNamedVerbose("Merge","Orient FB " << ffDist
                      << " " << fbDist
                      << " " << bfDist
                      << " " << bbDist);
        return kFrontBack;
    }
    else if (bfDist<fbDist && bfDist<ffDist && bfDist<bbDist) {
        // A back to front track.
        CaptNamedVerbose("Merge","Orient BF " << ffDist
                      << " " << fbDist
                      << " " << bfDist
                      << " " << bbDist);
        return kBackFront;
    }
    else if (bbDist<fbDist && bbDist<bfDist && bbDist<ffDist) {
        // A back to back track.
        CaptNamedVerbose("Merge","Orient BB " << ffDist
                      << " " << fbDist
                      << " " << bfDist
                      << " " << bbDist);
        return kBackBack;
    }

    // This can't happen!
    CaptError("This can't happen!");
    return kNotClose;
}

CP::THandle<CP::TReconTrack> 
CP::TMergeTracks::MergeTracks(CP::THandle<CP::TReconTrack> first, 
                              CP::THandle<CP::TReconTrack> second) {
    typedef std::deque< CP::THandle<CP::TReconCluster> > ClusterContainer;

    // Start by inserting the clusters from the first and second tracks into a
    // deque.  The second track could go in a vector, but the deque should be
    // OK.
    ClusterContainer clusters(first->GetNodes().size()
                              + second->GetNodes().size());
    clusters.clear();
    for (CP::TReconNodeContainer::iterator n = first->GetNodes().begin();
         n != first->GetNodes().end(); ++n) {
        CP::THandle<CP::TReconCluster> c = (*n)->GetObject();
        if (!c) {
            CaptError("Node not made from a cluster");
            throw;
        }
        clusters.push_back(c);
    }
            
    ClusterContainer secondClusters(second->GetNodes().size());
    secondClusters.clear();
    for (CP::TReconNodeContainer::iterator n = second->GetNodes().begin();
         n != second->GetNodes().end(); ++n) {
        CP::THandle<CP::TReconCluster> c = (*n)->GetObject();
        if (!c) {
            CaptError("Node not made from a cluster");
            throw;
        }
        // Check which clusters are overlapping.
        bool overlap = false;
        for (ClusterContainer::iterator s = clusters.begin();
             s != clusters.end(); ++s) {
            if (c == *s) {
                CaptDebug("Overlapping cluster at " 
                          << s-clusters.begin()
                          << " of " << clusters.size()-1);
                overlap = true;
                break;
            }
        }
        if (!overlap) secondClusters.push_back(c);
    }

    Orientation orient = TrackOrientation(first, second);

    // Figure out which ends are closest and use that to determine which order
    // the clusters get added.
    switch (orient) {
    case kFrontFront:
        // A front to front track.
        for (ClusterContainer::iterator s = secondClusters.begin();
             s != secondClusters.end(); ++s) {
            if (*s != clusters.front()) continue;
            clusters.push_front(*s);
        }
        break;
    case kFrontBack:
        // A front to back track.
        for (ClusterContainer::reverse_iterator s = secondClusters.rbegin();
             s != secondClusters.rend(); ++s) {
            if (*s != clusters.front()) continue;
            clusters.push_front(*s);
        }
        break;
    case kBackFront:
        // A back to front track.
        for (ClusterContainer::iterator s = secondClusters.begin();
             s != secondClusters.end(); ++s) {
            if (*s != clusters.back()) continue;
            clusters.push_back(*s);
        }
        break;
    case kBackBack:
        // A back to back track.
        for (ClusterContainer::reverse_iterator s = secondClusters.rbegin();
             s != secondClusters.rend(); ++s) {
            if (*s != clusters.back()) continue;
            clusters.push_back(*s);
        }
        break;
    default:
        CaptError("This can't happen!");
        return CP::THandle<CP::TReconTrack>();
    }

    // Create the new track.
    CP::THandle<CP::TReconTrack> track
        = CP::CreateTrack("TMergeTracks",clusters.begin(), clusters.end());

    TTrackFit fitter;
    track = fitter(track);

    return track;
}

double 
CP::TMergeTracks::MatchGoodness(CP::THandle<CP::TReconTrack> t1,
                                CP::THandle<CP::TReconTrack> t2) {
    Orientation orient = TrackOrientation(t1,t2);
    if (orient == kNotClose) return 1E+300;

    // The overall direction covariance.
    TMatrixD dirCov(3,3);

    // Get the matching position for the first track.  The direction is set so
    // that it points toward the "interior" of the track.
    TVector3 pos1;
    TMatrixD posCov1(3,3);
    TVector3 dir1;
    double q1 = 1.0;
    if (t1->GetNDOF() > 1) q1 = t1->GetQuality()/t1->GetNDOF();
    if (q1<1.0) q1 = 1.0;
    switch (orient) {
    case kFrontFront:
    case kFrontBack:
        pos1 = t1->GetFront()->GetPosition().Vect();
        dir1 = t1->GetFront()->GetDirection();
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                posCov1(i,j) = q1*t1->GetFront()->GetPositionCovariance(i,j);
                dirCov(i,j) += q1*t1->GetFront()->GetDirectionCovariance(i,j);
            }
        }
        break;
    case kBackFront:
    case kBackBack:
        pos1 = t1->GetBack()->GetPosition().Vect();
        dir1 = -t1->GetBack()->GetDirection();
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                posCov1(i,j) = q1*t1->GetBack()->GetPositionCovariance(i,j);
                dirCov(i,j) += q1*t1->GetBack()->GetDirectionCovariance(i,j);
            }
        }
        break;
    default:
        CaptError("This can't happen!");
    }
        

    // Get the matching position for the second track.  The direction is set
    // so that it points toward the "interior" of the track.
    TVector3 pos2;
    TMatrixD posCov2(3,3);
    TVector3 dir2;
    double q2 = 1.0;
    if (t2->GetNDOF() > 1) q2 = t2->GetQuality()/t2->GetNDOF();
    if (q2<1.0) q2 = 1.0;
    switch (orient) {
    case kFrontFront:
    case kBackFront:
        pos2 = t2->GetFront()->GetPosition().Vect();
        dir2 = t2->GetFront()->GetDirection();
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                posCov2(i,j) = q2*t2->GetFront()->GetPositionCovariance(i,j);
                dirCov(i,j) += q2*t2->GetFront()->GetDirectionCovariance(i,j);
            }
        }
        break;
    case kFrontBack:
    case kBackBack:
        pos2 = t2->GetBack()->GetPosition().Vect();
        dir2 = - t2->GetBack()->GetDirection();
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                posCov2(i,j) = q2*t2->GetBack()->GetPositionCovariance(i,j);
                dirCov(i,j) += q2*t2->GetBack()->GetDirectionCovariance(i,j);
            }
        }
        break;
    default:
        CaptError("This can't happen!");
    }

    // Make sure that the tracks don't overlap.  The directions are pointing
    // to the interior of the tracks, so they should be pointing in opposite
    // directions.
    double dirOverlap = dir1*dir2;

    // Make sure the two directions are pointing in the same general direction
    // (i.e. not opposites).
    double v = dir1*dir2;
    if (v<0) dir2 = -dir2;

    // Find out how well the directions match.
    TVector3 dirDiff = dir1-dir2;

    dirCov.InvertFast();
    double dirGoodness = dirDiff*(dirCov*dirDiff);

    // Find the average position between the two ends.
    TVector3 posAvg = 0.5*(pos1 + pos2);

    // Project the first track forward to the point of closest approach to the
    // average position and find out how well it matchs the average position.
    TVector3 pos1Diff = pos1-posAvg;
    double dist = pos1Diff*dir1;
    pos1Diff = pos1 - dist*dir1 - posAvg;
    posCov1.InvertFast();
    double pos1Goodness = pos1Diff*(posCov1*pos1Diff);

    // Find out how well the second track matches the average position.
    TVector3 pos2Diff = pos2-posAvg;
    dist = pos2Diff*dir2;
    pos2Diff = pos2 - dist*dir2 - posAvg;
    posCov2.InvertFast();
    double pos2Goodness = pos2Diff*(posCov2*pos2Diff);

    double result = dirGoodness + pos1Goodness + pos2Goodness;
    if (dirOverlap > 0.0) result += 1000.0;

    CaptNamedInfo("Merge","Goodness: " << result
                  << " Dir: " << dirGoodness
                  << " P1: " << pos1Goodness
                  << " P2: " << pos2Goodness
                  << " Ovlp: " << dirOverlap);

    return result;
}


CP::THandle<CP::TAlgorithmResult>
CP::TMergeTracks::Process(const CP::TAlgorithmResult& input,
                          const CP::TAlgorithmResult&,
                          const CP::TAlgorithmResult&) {
    
    CP::THandle<CP::TReconObjectContainer> inputObjects 
        = input.GetResultsContainer();

    CaptLog("TMergeTracks Process " << GetEvent().GetContext());

    if (!inputObjects) {
        CaptError("No input objects");
        return CP::THandle<CP::TAlgorithmResult>();
    }

    // Create the output containers.
    CP::THandle<CP::TAlgorithmResult> result = CreateResult();
    std::auto_ptr<CP::TReconObjectContainer> 
        final(new CP::TReconObjectContainer("final"));

    // Create a stack to keep tracks in.  When the stack is empty, all the
    // tracks that need to be merge have been merge.
    typedef std::list< CP::THandle<CP::TReconTrack> > TrackList;
    TrackList trackList;

    // Make a copy of all of the tracks in the input objects, and save the
    // non-tracks to the output.
    for (CP::TReconObjectContainer::iterator t = inputObjects->begin();
         t != inputObjects->end(); ++t) {
        CP::THandle<CP::TReconTrack> track = *t;
        if (!track) {
            final->push_back(*t);
            continue;
        }
        trackList.push_back(track);
    }

    // Pop a track off the stack and see if it should be merged.  If the track
    // is merged, the result is pushed back on the stack.  If it doesn't get
    // merge, the track get's pushed into the final object container.
    while (!trackList.empty()) {
        CP::THandle<CP::TReconTrack> track1 = trackList.front();
        trackList.pop_front();
        CaptNamedInfo("Merge",
                      "Track Stack: " << trackList.size()
                      << "    Track size: " << track1->GetNodes().size()
                      << "    UID: " << track1->GetUniqueID());

        // Don't use a very short track as a base for merging.  The tracks are
        // in order of number of nodes (decreasing), and a short track will
        // already be merged, or doesn't have a good partner.
        if (track1->GetNodes().size() < 3) {
            CaptNamedInfo("Merge", "Save short track  (" 
                          << track1->GetUniqueID() << ")");
            TTrackMassFit massFitter;
            CP::THandle<CP::TReconTrack> massTrack = massFitter(track1);
            if (massTrack) track1 = massTrack;
            final->push_back(track1);
            continue;
        }

        int index = 0;
        for (TrackList::iterator t = trackList.begin();
             t!=trackList.end(); ++t) {
            ++index;
            CP::THandle<CP::TReconTrack> track2 = *t;
            
            CaptNamedInfo("Merge", "Check Tracks"
                          << " w/ stack: " << trackList.size() 
                          << "  tracks: " << track1->GetUniqueID()
                          << ", " << track2->GetUniqueID()
                          << "  sizes: " 
                          << track1->GetNodes().size() 
                          << ", " << track2->GetNodes().size());

            double match = MatchGoodness(track1,track2);
            if (match > fGoodnessCut) continue;

            ///////////////////////////////////////////////////////
            // If we get here then the tracks should be merged.
            ///////////////////////////////////////////////////////

            // Remove the track iterator (the track is held in track2).
            trackList.erase(t);

            CP::THandle<CP::TReconTrack> merged = MergeTracks(track1,track2);

            // Merge the tracks.
            CaptNamedInfo("Merge", "Merge Tracks " << track1->GetUniqueID()
                          << ", " << track2->GetUniqueID()
                          << " into " << merged->GetUniqueID());

            // Put the new track back into the list.
            trackList.push_front(merged);

            // Clear out the track variable.
            track1 = CP::THandle<CP::TReconTrack>();

            // Don't continue the loop.
            break;
        }

        // If we get to the bottom of the loop with a track, then push it on
        // to the final objects.
        if (track1) {
            CaptNamedInfo("Merge", "Save Track");
            TTrackMassFit massFitter;
            CP::THandle<CP::TReconTrack> massTrack = massFitter(track1);
            if (massTrack) track1 = massTrack;
            final->push_back(track1);
        }

    }
    
    std::sort(final->begin(), final->end(), CompareReconObjects());

    result->AddResultsContainer(final.release());

    return result;
}
