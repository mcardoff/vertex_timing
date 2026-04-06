// export_training_data.cxx
//
// Exports one CSV row per cluster per passing event for DNN retraining.
//
// Columns
// -------
//   event_num            : global event index (chain read entry)
//   delta_z              : cluster_z - reco_vtx_z  [mm]  (precision-weighted centroid)
//   delta_z_resunits     : delta_z / cluster_z_sigma  (dimensionless)
//   cluster_z_sigma      : precision-weighted z uncertainty  [mm]
//   cluster_d0           : precision-weighted transverse impact parameter  [mm]
//   cluster_d0_sigma     : uncertainty on cluster_d0  [mm]
//   cluster_qOverP       : precision-weighted q/p  [MeV^{-1}]
//   cluster_qOverP_sigma : uncertainty on cluster_qOverP  [MeV^{-1}]
//   cluster_sumpt        : scalar sum of constituent track pT  [GeV]
//   cluster_time_sigma   : precision-weighted HGTD time uncertainty  [ps]  (= sigmas[0])
//   cluster_n_tracks     : number of constituent tracks (integer stored as float)
//   n_hs_tracks          : number of truth-matched HS (truthVtx_idx==0) tracks
//   delta_t              : cluster_time - recoVtxTime[0]  [ps]
//   label_purity         : 1 if cluster purity >= 0.5, else 0
//                          (select higher-purity clusters as signal)
//   purity               : raw HS pT fraction of the cluster  [0, 1]
//                          (continuous target for a regression / score network)
//
// Event selection mirrors clustering_dt / event_processing.h exactly:
//   passBasicCuts() + passJetPtCut(), track association at MAX_NSIGMA (3.0),
//   iterative clustering (distCut=DIST_CUT_CONE), real HGTD times,
//   checkValidTimes=true.
//
// The DNN is NOT run; features are computed raw so they can be used for
// retraining with any choice of normalisation.

#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>

#include <TROOT.h>
#include <TChain.h>
#include <TTreeReader.h>

#include "clustering_constants.h"
#include "clustering_structs.h"
#include "clustering_functions.h"
#include "event_processing.h"

using namespace MyUtl;

// ---------------------------------------------------------------------------
// Compute unnormalised features for a cluster (same arithmetic as
// calcFeatures() in clustering_structs.h, without normalisation so the
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
    float cluster_time_sigma;
    float cluster_n_tracks;
    float delta_t;
};

static RawFeatures calcRawFeatures(
    const Cluster& cluster,
    BranchPointerWrapper* branch
) {
    float znum  = 0.f, zden  = 0.f;
    float dnum  = 0.f, dden  = 0.f;
    float qpnum = 0.f, qpden = 0.f;
    float sumpt = 0.f;

    for (int trk : cluster.trackIndices) {
        float trk_z   = branch->trackZ0[trk], trk_var_z = branch->trackVarZ0[trk];
        float trk_d   = branch->trackD0[trk], trk_var_d = branch->trackVarD0[trk];
        float trk_q   = branch->trackQP[trk], trk_var_q = branch->trackVarQp[trk];
        float trk_pt  = branch->trackPt[trk];

        znum  += trk_z / trk_var_z;  zden  += 1.f / trk_var_z;
        dnum  += trk_d / trk_var_d;  dden  += 1.f / trk_var_d;
        qpnum += trk_q / trk_var_q;  qpden += 1.f / trk_var_q;
        sumpt += trk_pt;
    }

    float cluster_z       = znum  / zden;
    float cluster_z_sigma = 1.f   / std::sqrt(zden);
    float cluster_d0      = dnum  / dden;
    float cluster_d0_sig  = 1.f   / std::sqrt(dden);
    float cluster_qp      = qpnum / qpden;
    float cluster_qp_sig  = 1.f   / std::sqrt(qpden);
    float delta_z         = cluster_z - branch->recoVtxZ[0];
    float delta_z_ru      = delta_z / cluster_z_sigma;

    float delta_t = static_cast<float>(cluster.values.at(0)) - branch->recoVtxTime[0];

    return { delta_z, delta_z_ru, cluster_z_sigma,
             cluster_d0, cluster_d0_sig,
             cluster_qp, cluster_qp_sig,
             sumpt,
             static_cast<float>(cluster.sigmas.at(0)),
             static_cast<float>(cluster.trackIndices.size()),
             delta_t };
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

    csv << "event_num,"
        << "delta_z,delta_z_resunits,cluster_z_sigma,"
        << "cluster_d0,cluster_d0_sigma,"
        << "cluster_qOverP,cluster_qOverP_sigma,"
        << "cluster_sumpt,cluster_time_sigma,cluster_n_tracks,"
        << "delta_t,n_hs_tracks,"
        << "label_purity,purity\n";
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

        std::vector<int> tracks =
            getAssociatedTracks(&branch, MIN_TRACK_PT, MAX_TRACK_PT, MAX_NSIGMA);
        if (tracks.empty()) continue;

        // ---- Clustering: iterative, real HGTD times, calcPurityFlag=true ----
        std::vector<Cluster> clusters =
            clusterTracksInTime(
                tracks, &branch,
                DIST_CUT_CONE,
                /*useSmearedTimes*/ false,
                /*checkTimeValid*/  true,
                /*smearRes*/        IDEAL_TRACK_RES,
                ClusteringMethod::ITERATIVE,
                /*usez0*/           false,
                /*sortTracks*/      false,
                /*calcPurityFlag*/  true);

        if (clusters.empty()) continue;
        ++nPassEvent;

        // ---- Event-level: count forward HS tracks ----
        int nFTrack = 0, nFTrackHS = 0, nFTrackPU = 0;
        branch.countForwardTracks(nFTrack, nFTrackHS, nFTrackPU, tracks, /*checkTimeValid=*/true);

        // ---- One row per cluster ----
        for (const Cluster& cl : clusters) {
            if (cl.values.empty()) continue;
	    if (cl.nConstituents < 3) continue;
	    
            RawFeatures f    = calcRawFeatures(cl, &branch);
            float purity     = static_cast<float>(cl.purity);
            int label_purity = (purity >= 0.2f) ? 1 : 0;

            csv << evtNum                  << ","
                << f.delta_z               << ","
                << f.delta_z_resunits      << ","
                << f.cluster_z_sigma       << ","
                << f.cluster_d0            << ","
                << f.cluster_d0_sigma      << ","
                << f.cluster_qOverP        << ","
                << f.cluster_qOverP_sigma  << ","
                << f.cluster_sumpt         << ","
                << f.cluster_time_sigma    << ","
                << f.cluster_n_tracks      << ","
                << f.delta_t               << ","
                << nFTrackHS               << ","
                << label_purity            << ","
                << purity                  << "\n";
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
