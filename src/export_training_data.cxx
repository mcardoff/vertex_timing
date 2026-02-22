// export_training_data.cxx
//
// Exports one CSV row per cluster per passing event for DNN retraining.
//
// Columns
// -------
//   event_num            : global event index (chain read entry)
//   delta_t              : cluster_time - truth_HS_time  [ps]
//   delta_z              : cluster_z - reco_vtx_z  [mm]  (precision-weighted centroid)
//   delta_z_resunits     : delta_z / cluster_z_sigma  (dimensionless)
//   cluster_z_sigma      : precision-weighted z uncertainty  [mm]
//   cluster_d0           : precision-weighted transverse impact parameter  [mm]
//   cluster_d0_sigma     : uncertainty on cluster_d0  [mm]
//   cluster_qOverP       : precision-weighted q/p  [MeV^{-1}]
//   cluster_qOverP_sigma : uncertainty on cluster_qOverP  [MeV^{-1}]
//   cluster_sumpt        : scalar sum of constituent track pT  [GeV]
//   cluster_dR           : min deltaR between cluster eta/phi centroid and nearest topo jet
//   has_hs               : 1 if >=1 constituent track is truth-matched to the HS vertex
//                          (trackToTruthvtx[trk] == 0), else 0
//   n_hs_tracks          : number of forward HS tracks in the event with a valid HGTD time
//                          (event-level: identical for all clusters with the same event_num)
//   label                : 1 if the cluster passes the timing efficiency window
//                          (|delta_t| < 3*PASS_SIGMA = 60 ps) AND has_hs == 1, else 0
//
// Event selection mirrors clustering_dt / event_processing.h exactly:
//   passBasicCuts() + passJetPtCut(), track association at MAX_NSIGMA (3.0),
//   cone clustering (distCut=3.0, smearRes=10 ps), real HGTD times, checkValidTimes=true.
//
// The DNN is NOT run; features are computed raw so they can be used for retraining.

#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>

#include <TROOT.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TVector2.h>

#include "clustering_constants.h"
#include "clustering_structs.h"
#include "clustering_functions.h"
#include "event_processing.h"

using namespace MyUtl;

// ---------------------------------------------------------------------------
// Compute the pT-weighted eta/phi centroid of a cluster, then find the
// minimum deltaR to any reconstructed topo jet passing MIN_JETPT.
// ---------------------------------------------------------------------------
static double clusterMinDR(
    const Cluster& cluster,
    BranchPointerWrapper* branch
) {
    // Build pT-weighted eta/phi centroid from constituent tracks
    double sumPt = 0., etaSum = 0., phiSum = 0.;
    for (int trk : cluster.trackIndices) {
        double pt  = branch->trackPt[trk];
        double eta = branch->trackEta[trk];
        double phi = branch->trackPhi[trk];
        sumPt   += pt;
        etaSum  += pt * eta;
        phiSum  += pt * phi;
    }
    if (sumPt <= 0.) return -1.;
    double cEta = etaSum / sumPt;
    double cPhi = phiSum / sumPt;

    double minDR = 1e9;
    for (int j = 0; j < (int)branch->topoJetPt.GetSize(); ++j) {
        if (branch->topoJetPt[j] < MIN_JETPT) continue;
        double deta = branch->topoJetEta[j] - cEta;
        double dphi = TVector2::Phi_mpi_pi(branch->topoJetPhi[j] - cPhi);
        double dR   = std::sqrt(deta*deta + dphi*dphi);
        if (dR < minDR) minDR = dR;
    }
    return (minDR < 1e8) ? minDR : -1.;
}

// ---------------------------------------------------------------------------
// Compute the unnormalised features for a cluster (same arithmetic as
// calcFeatures() in clustering_structs.h, but without normalisation so the
// retraining script can choose its own scaler).
// ---------------------------------------------------------------------------
struct RawFeatures {
    float delta_z;
    float delta_z_resunits;
    float cluster_z_sigma;
    float cluster_d0;
    float cluster_d0_sigma;
    float cluster_qOverP;
    float cluster_qOverP_sigma;
    float cluster_sumpt;
};

static RawFeatures calcRawFeatures(
    const Cluster& cluster,
    BranchPointerWrapper* branch
) {
    float znum  = 0., zden  = 0.;
    float dnum  = 0., dden  = 0.;
    float qpnum = 0., qpden = 0.;
    float sumpt = 0.;

    for (int trk : cluster.trackIndices) {
        float trk_z     = branch->trackZ0[trk],   trk_var_z = branch->trackVarZ0[trk];
        float trk_d     = branch->trackD0[trk],   trk_var_d = branch->trackVarD0[trk];
        float trk_q     = branch->trackQP[trk],   trk_var_q = branch->trackVarQp[trk];
        float trk_pt    = branch->trackPt[trk];

        znum  += trk_z  / trk_var_z;  zden  += 1.f / trk_var_z;
        dnum  += trk_d  / trk_var_d;  dden  += 1.f / trk_var_d;
        qpnum += trk_q  / trk_var_q;  qpden += 1.f / trk_var_q;
        sumpt += trk_pt;
    }

    float cluster_z       = znum  / zden;
    float cluster_z_sigma = 1.f   / std::sqrt(zden);
    float cluster_d0      = dnum  / dden;
    float cluster_d0_sig  = 1.f   / std::sqrt(dden);
    float cluster_qp      = qpnum / qpden;
    float cluster_qp_sig  = 1.f   / std::sqrt(qpden);

    float ref_z      = branch->recoVtxZ[0];
    float delta_z    = cluster_z - ref_z;
    float delta_z_ru = delta_z   / cluster_z_sigma;

    return { delta_z, delta_z_ru, cluster_z_sigma,
             cluster_d0, cluster_d0_sig,
             cluster_qp, cluster_qp_sig,
             sumpt };
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
auto main() -> int {
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;

    // --- Data source (same as clustering_dt) ---
    TChain chain("ntuple");
    setupChain(chain, "../../ntuple-hgtd/");
    TTreeReader reader(&chain);
    BranchPointerWrapper branch(reader);

    // --- Output CSV ---
    const char* outPath = "../cluster_data_for_training.csv";
    std::ofstream csv(outPath);
    if (!csv.is_open()) {
        std::cerr << "ERROR: cannot open " << outPath << " for writing\n";
        return 1;
    }

    // Write header
    csv << "event_num,delta_t,delta_z,delta_z_resunits,"
        << "cluster_z_sigma,cluster_d0,cluster_d0_sigma,"
        << "cluster_qOverP,cluster_qOverP_sigma,cluster_sumpt,"
        << "cluster_dR,has_hs,n_hs_tracks,label\n";
    csv << std::setprecision(6) << std::fixed;

    // --- Event loop ---
    const Long64_t nTotal = chain.GetEntries();
    Long64_t nPassEvent   = 0;
    Long64_t nRowsWritten = 0;

    std::cout << "Starting export loop over " << nTotal << " events\n";

    while (reader.Next()) {
        const Long64_t evtNum = chain.GetReadEntry();
        if ((evtNum + 1) % 5000 == 0)
            std::cout << "Progress: " << evtNum + 1 << " / " << nTotal << "\n";

        // ---- Event selection (mirrors processEventData) ----
        if (!branch.passBasicCuts()) continue;
        if (!branch.passJetPtCut())  continue;

        // Track association at 3.0σ (MAX_NSIGMA, same as clustering_dt after
        // the linter change restoring MAX_NSIGMA = 3.0)
        std::vector<int> tracks =
            getAssociatedTracks(&branch, MIN_TRACK_PT, MAX_TRACK_PT, MAX_NSIGMA);

        if (tracks.empty()) continue;

        // Count forward HS tracks in the event (same 3σ list, checkValidTimes=true).
        // This is an event-level quantity — all clusters in the same event share it.
        int nFTrack = 0, nFTrackHS = 0, nFTrackPU = 0;
        branch.countForwardTracks(nFTrack, nFTrackHS, nFTrackPU, tracks, /*checkTimeValid=*/true);

        // ---- Clustering: real HGTD times, checkValidTimes=true, cone algorithm ----
        // calcPurityFlag=true so purity is available for has_hs; no ML scoring needed.
        std::vector<Cluster> clusters =
            clusterTracksInTime(
                tracks, &branch,
                /* distanceCut */     3.0,
                /* useSmearedTimes */ false,
                /* checkTimeValid */  true,
                /* smearRes */        10.0,
                /* useCone */         true,
                /* usez0 */           false,
                /* calcPurityFlag */  true
            );

        if (clusters.empty()) continue;
        ++nPassEvent;

        const float truthTime = branch.truthVtxTime[0];

        // ---- One row per cluster ----
        for (const Cluster& cl : clusters) {
            if (cl.values.empty()) continue;

            // delta_t: cluster time - truth HS time
            float cluster_time = static_cast<float>(cl.values[0]);
            float delta_t      = cluster_time - truthTime;

            // Unnormalised kinematic features
            RawFeatures f = calcRawFeatures(cl, &branch);

            // Min dR to nearest topo jet
            double dR = clusterMinDR(cl, &branch);

            // has_hs: any constituent track truth-matched to HS vertex (idx == 0)
            int has_hs = 0;
            for (int trk : cl.trackIndices) {
                if (branch.trackToTruthvtx[trk] == 0) { has_hs = 1; break; }
            }

            // label: cluster time is within timing efficiency window of truth HS
            int label = (has_hs && (std::abs(delta_t) < 3.f * static_cast<float>(PASS_SIGMA))) ? 1 : 0;

            csv << evtNum              << ","
                << delta_t             << ","
                << f.delta_z           << ","
                << f.delta_z_resunits  << ","
                << f.cluster_z_sigma   << ","
                << f.cluster_d0        << ","
                << f.cluster_d0_sigma  << ","
                << f.cluster_qOverP    << ","
                << f.cluster_qOverP_sigma << ","
                << f.cluster_sumpt     << ","
                << dR                  << ","
                << has_hs              << ","
                << nFTrackHS           << ","
                << label               << "\n";
            ++nRowsWritten;
        }
    }

    csv.close();
    std::cout << "\nDone.\n"
              << "  Passing events : " << nPassEvent   << " / " << nTotal << "\n"
              << "  Rows written   : " << nRowsWritten << "\n"
              << "  Output file    : " << outPath      << "\n";

    return 0;
}
