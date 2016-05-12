#include "TSplitTracks.hxx"
#include "TTrackFit.hxx"
#include "ClusterDistance.hxx"
#include "CreateTrack.hxx"
#include "CompareReconObjects.hxx"

#include <THandle.hxx>
#include <TReconTrack.hxx>
#include <TReconCluster.hxx>
#include <TCaptLog.hxx>
#include <HEPUnits.hxx>
#include <TUnitsTable.hxx>
#include <TRuntimeParameters.hxx>
#include <ostreamTVector3.hxx>

#include <TMatrixD.h>
#include <TPrincipal.h>

#include <memory>
#include <cmath>
#include <sstream>

CP::TSplitTracks::TSplitTracks()
    : TAlgorithm("TSplitTracks", 
                 "Split Tracks with Kinks and Gaps") {\
    fThreeInLineCut = 3.0;
    fRadiusOfCurvature = 25*unit::cm;
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
        CaptNamedInfo("Split", "Save track clusters.  Too short");
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
        = CP::CreateTrackFromClusters("TSplitTracks",begin, end);

    if (!track) {
        CaptNamedInfo("Split", "Track not created by CreateTrackFromClusters");
    }
    
    CP::TTrackFit fitter;
    track = fitter(track);
    
    CaptNamedInfo("Split", "Save track with "
                  << end-begin << " clusters"
                  << " (UID " << track->GetUniqueID() << ")");
    
    container.push_back(track);
}

double
CP::TSplitTracks::RadiusOfCurvature(CP::THandle<CP::TReconCluster> a,
                                    CP::THandle<CP::TReconCluster> b, 
                                    CP::THandle<CP::TReconCluster> c) {
    // The assumption is that B lies between A and C, but the code will work
    
    double len = (c->GetPosition().Vect()-a->GetPosition().Vect()).Mag();
    TVector3 v1 = (b->GetPosition().Vect()-a->GetPosition().Vect()).Unit();
    TVector3 v2 = (c->GetPosition().Vect()-b->GetPosition().Vect()).Unit();
    double csgn = v1*v2;
    // If this is a kink of more than 90 deg, then set the radius to zero.
    if (csgn < 0.0) return 0.0;
    double sgn = std::sqrt(1.0 - csgn*csgn);
    return len / (2.0*sgn);
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

    // Add in a factor for multiple scattering
    double radLen = 14*unit::cm; // For liquid argon.
    double X = std::abs(proj)/radLen;
    double P = 100*unit::MeV; // HACK! 
    if (X < 0.001) X = 0.001;
    // Set the minimum amount of scattering.  This isn't very physical, but
    // there needs to be scattering or the fit doesn't work.
    double scatter = (1.0+0.038*std::log(X))*std::sqrt(X)*(13.6*unit::MeV)/(P);
    scatter = std::abs(proj*scatter);
    for (int i=0; i<3; ++i) {
        cov(i,i) += scatter;
        cov(i,i) += scatter;
        cov(i,i) += scatter;
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

    double p0[3] = {0.0,0.0,0.0};
    double p1[3] = {1.0,0.0,0.0};
    double p2[3] = {0.0,1.0,0.0};
    double p3[3] = {0.0,0.0,1.0};
    double x[3];

    // Find out where to start the fit in the backward direction.  This makes
    // sure the distance isn't too big.
    ClusterContainer::iterator startStep = here;
    while (startStep != begin) {
        --startStep;
        if (here-startStep+1 < minStep) continue;
        double r = ((*startStep)->GetPosition().Vect()
                    -(*here)->GetPosition().Vect()).Mag();
        if (r > 4*unit::cm) break;
    }

    // Do a simple fit in the backward direction.  This includes the cluster
    // being checked.
    TPrincipal pca1(3,"");
    TPrincipal pca1e(3,"");
    for (ClusterContainer::iterator i = startStep; i != here+1; ++i) {
        double row[3] = {(*i)->GetPosition().X(),
                         (*i)->GetPosition().Y(),
                         (*i)->GetPosition().Z()};
        pca1.AddRow(row);
        pca1e.AddRow(row);
        TVector3 corner;

        corner = (*i)->GetPosition().Vect() + 1.5*(*i)->GetLongAxis();
        double row1[3] = {corner.X(), corner.Y(), corner.Z()};
        pca1e.AddRow(row1);

        corner = (*i)->GetPosition().Vect() - 1.5*(*i)->GetLongAxis();
        double row2[3] = {corner.X(), corner.Y(), corner.Z()};
        pca1e.AddRow(row2);

        corner = (*i)->GetPosition().Vect() + (*i)->GetMajorAxis();
        double row3[3] = {corner.X(), corner.Y(), corner.Z()};
        pca1e.AddRow(row3);

        corner = (*i)->GetPosition().Vect() - (*i)->GetMajorAxis();
        double row4[3] = {corner.X(), corner.Y(), corner.Z()};
        pca1e.AddRow(row4);

        corner = (*i)->GetPosition().Vect() + 1.5*(*i)->GetMinorAxis();
        double row5[3] = {corner.X(), corner.Y(), corner.Z()};
        pca1e.AddRow(row5);

        corner = (*i)->GetPosition().Vect() - 1.5*(*i)->GetMinorAxis();
        double row6[3] = {corner.X(), corner.Y(), corner.Z()};
        pca1e.AddRow(row6);
    }
    pca1.MakePrincipals();
    pca1e.MakePrincipals();
    pca1.P2X(p0,x,3);
    TVector3 base1(x);
    pca1.P2X(p1,x,3);
    TVector3 dir1(x); 
    dir1 = (dir1-base1).Unit();
    TVector3 dirSense1
        = (*here)->GetPosition().Vect() - (*startStep)->GetPosition().Vect();
    double sense1 = dir1*dirSense1;
    if (sense1 < 0.0) dir1 = - dir1;
    
    // Find out where to end the fit in the forward direction.
    ClusterContainer::iterator endStep = here;
    while (endStep != end) {
        if (end-here < minStep) continue;
        double r = ((*endStep)->GetPosition().Vect()
                    -(*here)->GetPosition().Vect()).Mag();
        ++endStep;
        if (r > 4*unit::cm) break; // Yes... the "if" comes after the increment.
    }

    // Do a simple fit in the forward direction.  This includes the cluster
    // being checked.
    int foreStep = end-here;
    if (foreStep > maxStep) foreStep = maxStep;
    TPrincipal pca2(3,"");
    TPrincipal pca2e(3,"");
    for (ClusterContainer::iterator i = here; i != endStep; ++i) {
        double row[3] 
            = {(*i)->GetPosition().X(),
               (*i)->GetPosition().Y(),
               (*i)->GetPosition().Z()};
        pca2.AddRow(row);
        pca2e.AddRow(row);
        TVector3 corner;

        corner = (*i)->GetPosition().Vect() + 1.5*(*i)->GetLongAxis();
        double row1[3] = {corner.X(), corner.Y(), corner.Z()};
        pca2e.AddRow(row1);

        corner = (*i)->GetPosition().Vect() - 1.5*(*i)->GetLongAxis();
        double row2[3] = {corner.X(), corner.Y(), corner.Z()};
        pca2e.AddRow(row2);

        corner = (*i)->GetPosition().Vect() + 1.5*(*i)->GetMajorAxis();
        double row3[3] = {corner.X(), corner.Y(), corner.Z()};
        pca2e.AddRow(row3);

        corner = (*i)->GetPosition().Vect() - 1.5*(*i)->GetMajorAxis();
        double row4[3] = {corner.X(), corner.Y(), corner.Z()};
        pca2e.AddRow(row4);

        corner = (*i)->GetPosition().Vect() + 1.5*(*i)->GetMinorAxis();
        double row5[3] = {corner.X(), corner.Y(), corner.Z()};
        pca2e.AddRow(row5);

        corner = (*i)->GetPosition().Vect() - 1.5*(*i)->GetMinorAxis();
        double row6[3] = {corner.X(), corner.Y(), corner.Z()};
        pca2e.AddRow(row6);
    }
    pca2.MakePrincipals();
    pca2e.MakePrincipals();
    pca2.P2X(p0,x,3);
    TVector3 base2(x);
    pca2.P2X(p1,x,3);
    TVector3 dir2(x); 
    dir2 = (dir2-base2).Unit();
    TVector3 dirSense2
        = (*(endStep-1))->GetPosition().Vect() - (*here)->GetPosition().Vect();
    double sense2 = dir2*dirSense2;
    if (sense2 < 0.0) dir2 = - dir2;

    // Find the cosine of the angle between the fore and back segments.
    double dcos = dir1*dir2;
    dcos = std::max(-1.0,std::min(dcos,1.0));
    double angle = std::acos(dcos);

    // Find the error on the angle.  This changes the direction to make the
    // calculation.
    TVector3 dirDiff = dir2;
    if (dcos < 0) dirDiff = -dir2;
    dirDiff = (dirDiff - dir1).Unit();
    
    double transMag = 0.0;
    {
        pca1e.P2X(p0,x,3);
        TVector3 b(x[0],x[1],x[2]);
        pca1e.P2X(p1,x,3);
        TVector3 v(x[0],x[1],x[2]);
        v = (*pca1e.GetSigmas())(0)*(v - b);
        transMag += std::abs(v*dirDiff);
    }
    {
        pca1e.P2X(p0,x,3);
        TVector3 b(x[0],x[1],x[2]);
        pca1e.P2X(p2,x,3);
        TVector3 v(x[0],x[1],x[2]);
        v = (*pca1e.GetSigmas())(1)*(v - b);
        transMag += std::abs(v*dirDiff);
    }
    {
        pca1e.P2X(p0,x,3);
        TVector3 b(x[0],x[1],x[2]);
        pca1e.P2X(p3,x,3);
        TVector3 v(x[0],x[1],x[2]);
        v = (*pca1e.GetSigmas())(2)*(v - b);
        transMag += std::abs(v*dirDiff);
    }
    double dAng1 = atan2(transMag,(base1-(*here)->GetPosition().Vect()).Mag());
    
    transMag = 0.0;
    {
        pca2e.P2X(p0,x,3);
        TVector3 b(x[0],x[1],x[2]);
        pca2e.P2X(p1,x,3);
        TVector3 v(x[0],x[1],x[2]);
        v = (*pca2e.GetSigmas())(0)*(v - b);
        transMag += std::abs(v*dirDiff);
    }
    {
        pca2e.P2X(p0,x,3);
        TVector3 b(x[0],x[1],x[2]);
        pca2e.P2X(p2,x,3);
        TVector3 v(x[0],x[1],x[2]);
        v = (*pca2e.GetSigmas())(1)*(v - b);
        transMag += std::abs(v*dirDiff);
    }
    {
        pca2e.P2X(p0,x,3);
        TVector3 b(x[0],x[1],x[2]);
        pca2e.P2X(p3,x,3);
        TVector3 v(x[0],x[1],x[2]);
        v = (*pca2e.GetSigmas())(2)*(v - b);
        transMag += std::abs(v*dirDiff);
    }
    double dAng2 = atan2(transMag,(base2-(*here)->GetPosition().Vect()).Mag());

    double angle1 = angle - std::sqrt(dAng1*dAng1 + dAng2*dAng2);
    if (angle1 < 0.0) angle1 = 0.0;
    
    std::cout << "ANGLE " << angle
              << " " << angle1
              << " " << dAng1
              << " " << dAng2
              << std::endl;

    // Return the angle between the segments.
    return angle1;
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

    // Create the result container.
    CP::THandle<CP::TAlgorithmResult> result = CreateResult();

    // Create an output queue to stage output objects into.  This might
    // contain duplicate objects so it's filtered just before saving the
    // objects to final.
    CP::TReconObjectContainer outputStack;

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
            CaptNamedInfo("Split", "Non-track object in input"
                          << " UID: " << (*t)->GetUniqueID());
            outputStack.push_back(*t);
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
                      << "    Track size: " << track->GetNodes().size()
                      << "    Pop UID: " << track->GetUniqueID());
        
        // Make a copy of the clusters in the track.  These are basically the
        // track nodes.
        ClusterContainer clusters;
        std::ostringstream clusterIds;
        for (CP::TReconNodeContainer::iterator n = track->GetNodes().begin();
             n != track->GetNodes().end(); ++n) {
            CP::THandle<CP::TReconCluster> c = (*n)->GetObject();
            if (!c) {
                CaptError("Track node not made from a cluster");
                throw;
            }
            if (n != track->GetNodes().begin()) clusterIds << ", ";
            clusterIds << c->GetUniqueID();
            clusters.push_back(c);
        }
        CaptNamedInfo("Split","Clusters: " << clusterIds.str());
        
        ClusterContainer::iterator begin, end;
        begin = clusters.begin();
        end = clusters.end();

        // Protect against small tracks...
        if (end-begin < 3) {
            SaveTrack(outputStack,begin,end);
            continue;
        }

        // Check to see if the front clusters should be removed from the
        // track.  Add any clusters that are removed to the final objects.
        while ((begin + 3) != end) {
            double v = ThreeInLine(*begin, *(begin+1), *(begin+2));
            double r = RadiusOfCurvature(*begin, *(begin+1), *(begin+2));
            double d = ClusterDistance(**begin, **(begin+1));
            CaptNamedInfo("Split",
                          "Check front UID "
                          << (*begin)->GetUniqueID() << " --"
                          << " In Line: " << v
                          << " Radius: " << r
                          << " End Dist: " << d
                          << " Clusters: " << (end-begin)-1);
            if (v < fThreeInLineCut && d < fEndDistanceCut) break;
            if (r > fRadiusOfCurvature && d < fEndDistanceCut) break;
            CaptNamedInfo("Split",
                          "Discard front UID "
                          << (*begin)->GetUniqueID() << " --"
                          << " In Line: " << v
                          << " Radius: " << r
                          << " End Dist: " << d
                          << " Clusters: " << (end-begin)-1);
            outputStack.push_back(*begin);
            ++begin;
        }

        // Check to see if the back clusters should be removed from the track.
        // Add any clusters that are removed to the final objects.
        while (end != begin+3) {
            double v = ThreeInLine(*(end-1), *(end-2), *(end-3));
            double r = RadiusOfCurvature(*(end-1), *(end-2), *(end-3));
            double d = ClusterDistance(**(end-1), **(end-2));
            CaptNamedInfo("Split",
                          "Check back UID " 
                          << (*(end-1))->GetUniqueID() << " --"
                          << " In Line: " << v
                          << " Radius: " << r
                          << " End Dist: " << d
                          << " Clusters: " << (end-begin)-1);
            if (v < fThreeInLineCut && d < fEndDistanceCut) break;
            if (r > fRadiusOfCurvature && d < fEndDistanceCut) break;
            CaptNamedInfo("Split",
                          "Discard back UID " 
                          << (*(end-1))->GetUniqueID() << " --"
                          << " In Line: " << v
                          << " Radius: " << r
                          << " End Dist: " << d
                          << " Clusters: " << (end-begin)-1);
            outputStack.push_back(*(end-1));
            --end;
        }

        // Check to see if the track length has been reduced so much that it
        // should be split up into it's clusters.  The track will either be
        // added to the final set of objects, or split into clusters which are
        // added to the final set of objects.
        if (end-begin < 5) {
            if (end-begin < 3) {
                // This is a very short track stub after removing clusters, so
                // just save the clusters.
                std::copy(begin,end,std::back_inserter(outputStack));
                CaptNamedInfo("Split","Discard short track.  Length: "
                        << end-begin);
            }
            else if (end-begin == 3) {
                // These are three clusters in a row, so check if they are
                // consistent with a line.  If they are, make a track from the
                // clusters and save it to the final objects.
                double v = ThreeInLine(*begin, *(begin+1), *(begin+2));
                double r = RadiusOfCurvature(*begin, *(begin+1), *(begin+2));
                double d = std::max(ClusterDistance(**begin, **(begin+1)),
                                    ClusterDistance(**(begin+1), **(begin+2)));
                if (v < fThreeInLineCut && d < fEndDistanceCut) {
                    SaveTrack(outputStack,begin,end);
                } 
                else if (r > fRadiusOfCurvature && d < fEndDistanceCut) {
                    SaveTrack(outputStack,begin,end);
                } 
                else {
                    std::copy(begin,end,std::back_inserter(outputStack));
                    CaptNamedInfo("Split","Discard a short track."
                                  << " In Line: " << v
                                  << " Radius: " << r
                                  << " Gap : " << d
                                  << " Clusters " << end-begin);
                }
            }
            else if (end-begin == 4) {
                // These are four clusters in a row, so check if they are
                // consistent with a line.  If they are, make a track from the
                // clusters and save it to the final objects.
                double v1 = ThreeInLine(*begin, *(begin+1), *(begin+2));
                double v2 = ThreeInLine(*(begin+1), *(begin+2), *(begin+3));
                double v = std::max(v1,v2);
                double r1 = RadiusOfCurvature(*begin, *(begin+1), *(begin+2));
                double r2 = RadiusOfCurvature(*(begin+1),*(begin+2),*(begin+3));
                double r = std::min(r1, r2);
                double d = std::max(ClusterDistance(**begin, **(begin+1)),
                                    ClusterDistance(**(begin+2), **(begin+3)));
                if (v < fThreeInLineCut && d < fEndDistanceCut) {
                    SaveTrack(outputStack,begin,end);
                } 
                else if (r > fRadiusOfCurvature && d < fEndDistanceCut) {
                    SaveTrack(outputStack,begin,end);
                } 
                else {
                    std::copy(begin,end,std::back_inserter(outputStack));
                    CaptNamedInfo("Split","Discard a short track."
                                  << " In Line1: " << v1
                                  << " In Line2: " << v2
                                  << " Radius1: " << r1
                                  << " Radius2: " << r2
                                  << " Gap : " << d
                                  << " Clusters " << end-begin);
                }
            }

            // Track is disposed of so don't check it any more.
            continue;
        }

        CaptNamedInfo("Split", "Check long track"
                      << " (UID " << track->GetUniqueID() << ")");

        // Look for parts of the track where three clusters are not in a line.
        // The segment with the worst chi2 will be a "short" kink in the
        // track.  There are at least five clusters in the track when this bit
        // of code starts.  The kink can't be right at either end.
        ClusterContainer::iterator sharpestKink = begin+1;
        double worstChi2 = 0.0;
        double worstRadius = 1000.0*unit::mm;
        for (ClusterContainer::iterator i = begin+1; i+3 != end; ++i) {
            double r = RadiusOfCurvature(*(i),*(i+1), *(i+2));
            if (r > fRadiusOfCurvature) continue;
            double v = ThreeInLine(*(i),*(i+1), *(i+2));
            if (worstChi2 < v) {
                worstChi2 = v;
                worstRadius = r;
                sharpestKink = i+1;
            }
        }

        CaptNamedInfo("Split", "Sharpest kink found has chi2 " << worstChi2
                      << " (" << fThreeInLineCut << ")"
                      << " radius " << worstRadius
                      << " at " << sharpestKink-begin
                      << " of " << end-begin
                      << "  UID: " << sharpestKink->GetUniqueID());

        // If the worstChi2 is too big, there is a kink and the track should
        // be split.  The two pieces are put back on the stack and we start
        // over again.
        if (worstChi2 > fThreeInLineCut) {
            track = CreateTrackFromClusters("TSplitTracks",begin,sharpestKink+1);
            CaptNamedInfo("Split", "Stack track "
                          << " (UID " << track->GetUniqueID() << ")");
            trackStack.push_back(track);
            track = CreateTrackFromClusters("TSplitTracks",sharpestKink,end);
            CaptNamedInfo("Split", "Stack track "
                          << " (UID " << track->GetUniqueID() << ")");
            trackStack.push_back(track);
            CaptNamedInfo("Split","Split in line --"
                         << " X2: " << worstChi2
                         << " R: " << worstRadius
                         << "  Original: "   << end-begin
                         << "  First: "   << sharpestKink-begin+1
                         << "  Second: "   << end-sharpestKink);
            continue;
        }

        // Find the biggest gap in the track.  The can't be right at the end
        // of the track.  There are at least 5 clusters in the track to get to
        // this point, so the check (i+3 != end) should be OK.
        ClusterContainer::iterator biggestGap = begin+1;
        double maxDist = 0.0;
        for (ClusterContainer::iterator i = begin+1; i+3 != end; ++i) {
            double v = ClusterDistance(**(i),**(i+1));
            if (maxDist < v) {
                maxDist = v;
                biggestGap = i+1;
            }
        }

        CaptNamedInfo("Split", "Biggest gap " << maxDist
                      << " (" << fSplitDistanceCut << ")"
                      << " at " << biggestGap-begin
                      << " of " << end-begin
                      << "  UID: " << biggestGap->GetUniqueID());

        // If the biggest gap is too big the track should be split.  The two
        // pieces are put back on the stack and we start over again.
        if (maxDist > fSplitDistanceCut) {
            track = CreateTrackFromClusters("TSplitTracks",begin,biggestGap);
            CaptNamedInfo("Split", "Stack track with "
                          << " (UID " << track->GetUniqueID() << ")");
            trackStack.push_back(track);
            track = CreateTrackFromClusters("TSplitTracks",biggestGap,end);
            CaptNamedInfo("Split", "Stack track with "
                          << " (UID " << track->GetUniqueID() << ")");
            trackStack.push_back(track);
            CaptNamedInfo("Split","Split at gap --"
                         << "  Gap: " << maxDist 
                         << "  Original: "   << end-begin
                         << "  First: "   << biggestGap-begin
                         << "  Second: "   << end-biggestGap);
            continue;
        }

        // Look through the clusters and find the one with the sharpest kink
        // based on the track segments to either side of the current cluster.
        // This reuses the sharpestKink iterator.
        double biggestAngle = 0.0;
        ClusterContainer::iterator biggestBend = begin;
        for (ClusterContainer::iterator i = begin; i != end; ++i) {
            double v = KinkAngle(i,begin,end);
            if (biggestAngle < v) {
                biggestAngle = v;
                biggestBend = i;
            }
            CaptNamedInfo("SplitAngle","  Angle at " << i-begin
                          << " is " << v/unit::degree << " deg"
                          << " biggest angle is " << biggestAngle/unit::degree
                          << " at " << biggestBend-begin);
        }

        CaptNamedInfo("Split", "Biggest bend " << biggestAngle
                      << " (" << fKinkAngleCut << ")"
                      << " at " << biggestBend-begin
                      << " of " << end-begin
                      << "  UID: " << biggestBend->GetUniqueID());

        // If the kink angle is to big there is a kink and the track should
        // be split.  The two pieces are put back on the stack and we start
        // over again.
        if (biggestAngle > fKinkAngleCut) {
            track = CreateTrackFromClusters("TSplitTracks",begin,biggestBend+1);
            CaptNamedInfo("Split", "Stack track "
                          << " (UID " << track->GetUniqueID() << ")");
            trackStack.push_back(track);
            track = CreateTrackFromClusters("TSplitTracks",biggestBend,end);
            CaptNamedInfo("Split", "Stack track "
                          << " (UID " << track->GetUniqueID() << ")");
            trackStack.push_back(track);
            CaptNamedInfo("Split","Split at kink --"
                         << " Angle: " << biggestAngle
                         << "  Original: "   << end-begin
                         << "  First: "   << biggestBend-begin+1
                         << "  Second: "   << end-biggestBend);
            continue;
        }

        // Whatever gets here should just be saved...
        SaveTrack(outputStack,begin,end);
    }
    
    // Copy the outputStack objects to final.
    std::unique_ptr<CP::TReconObjectContainer> 
        final(new CP::TReconObjectContainer("final"));

    // A stack of clusters that have been seen.
    CP::TReconObjectContainer clusterStack;

    // Filter the tracks out of the outputStack and save their clusters to
    // the clusterStack.
    for (CP::TReconObjectContainer::iterator o = outputStack.begin();
         o != outputStack.end(); ++o) {
        CP::THandle<CP::TReconTrack> track = *o;
        if (track) {
            final->push_back(track);
            // Add the clusters in the track to the cluster stack.
            for (CP::TReconNodeContainer::iterator n 
                     = track->GetNodes().begin();
                 n != track->GetNodes().end(); ++n) {
                CP::THandle<CP::TReconCluster> c = (*n)->GetObject();
                if (!c) {
                    CaptError("Track node not made from a cluster");
                    throw;
                }
                clusterStack.push_back(c);
            }
            continue;
        }
        // Put any non cluster objects into the final container.
        CP::THandle<CP::TReconCluster> cluster = *o;
        if (cluster) continue;
        final->push_back(*o);
    }

    // Filter the clusters out of the outputStack and save any that are not in
    // the clusterStack to the final output container.
    for (CP::TReconObjectContainer::iterator o = outputStack.begin();
         o != outputStack.end(); ++o) {
        CP::THandle<CP::TReconCluster> cluster = *o;
        if (!cluster) continue;
        CP::TReconObjectContainer::iterator obj 
            = std::find(clusterStack.begin(), clusterStack.end(), cluster);
        if (obj != clusterStack.end()) continue;
        final->push_back(cluster);
        clusterStack.push_back(cluster);
    }

    std::sort(final->begin(), final->end(), CompareReconObjects());
    result->AddResultsContainer(final.release());

    return result;
}
