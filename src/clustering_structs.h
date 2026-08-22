#ifndef STRUCTS_H
#define STRUCTS_H

// ---------------------------------------------------------------------------
// clustering_structs.h
//   Core data structures shared across the entire analysis pipeline.
//   Two structs are defined here:
//
//   BranchPointerWrapper
//     Owns all TTreeReaderArray / TTreeReaderValue handles for a single
//     TTreeReader session.  Also hosts all per-event selection predicates
//     and kinematic helper methods so that event-loop code in
//     event_processing.h stays concise.
//
//   Cluster
//     Represents one time-space cluster produced by the clustering
//     algorithms in clustering_functions.h.  Stores the weighted-mean
//     time (and optionally z₀), per-constituent track indices, all raw
//     times (for spread diagnostics), and a score map keyed on Score.
//     Also hosts the ML feature extraction and score-update logic so
//     that the cluster is self-contained.
// ---------------------------------------------------------------------------

#include "clustering_includes.h"
#include "clustering_constants.h"
#include <cstdlib>

namespace MyUtl {

  // ---------------------------------------------------------------------------
  // BranchPointerWrapper
  //   Aggregates every TTreeReaderArray / TTreeReaderValue needed by the
  //   analysis into a single object that is passed by pointer throughout the
  //   event loop.  The constructor binds each member to its branch name so
  //   that branch-name strings appear only once in the codebase.
  //
  //   Selection methods are grouped as follows:
  //     Basic event cuts         — passBasicCuts, passJetPtCut
  //     Track / jet counters     — countForwardJets, countForwardTracks
  //     Testing-only HS filter   — passForwardHsTracks
  //     Track scoring helpers    — calcJetptDRScore, calcTrkptDRScore
  //     VBF signal region cuts   — see dedicated sub-section below
  // ---------------------------------------------------------------------------
  struct BranchPointerWrapper {
    TTreeReader& reader;

    TTreeReaderValue<float> weight;

    TTreeReaderArray<float> trackZ0;
    TTreeReaderArray<float> trackD0;
    TTreeReaderArray<float> trackPt;
    TTreeReaderArray<float> trackQP;
    TTreeReaderArray<float> trackEta;
    TTreeReaderArray<float> trackPhi;
    TTreeReaderArray<float> trackTheta;
    TTreeReaderArray<float> trackVarZ0;
    TTreeReaderArray<float> trackVarD0;
    TTreeReaderArray<float> trackVarQp;
    TTreeReaderArray<float> trackVarTheta;
    TTreeReaderArray<float> trackVarPhi;
    TTreeReaderArray<float> trackTime;
    TTreeReaderArray<float> trackTimeRes;
    TTreeReaderArray<int>   trackTimeValid;
    TTreeReaderArray<int>   trackToTruthvtx;
    TTreeReaderArray<int>   trackToParticle;
    TTreeReaderArray<bool>  trackQuality;
    TTreeReaderArray<int>   trackHgtdHits;
    TTreeReaderArray<int>   trackPrimHits;

    TTreeReaderArray<float> truthVtxZ;
    TTreeReaderArray<float> truthVtxTime;
    TTreeReaderArray<bool>  truthVtxIshs;

    TTreeReaderArray<float> recoVtxZ;
    TTreeReaderArray<float> recoVtxTime;
    TTreeReaderArray<float> recoVtxTimeRes;
    TTreeReaderArray<int>   recoVtxValid;

    TTreeReaderArray<float> topoJetPt;
    TTreeReaderArray<float> topoJetEta;
    TTreeReaderArray<float> topoJetPhi;

    TTreeReaderArray<float> truthHSJetPt;
    TTreeReaderArray<float> truthHSJetEta;
    TTreeReaderArray<float> truthHSJetPhi;

    // NOTE: TruthITPUJet_* / TruthOOTPUJet_* were bound here but never read by
    // any active code path, so the bindings are gone. Binding is not free: a
    // TTreeReaderArray on a branch the tree lacks is an error even if the value
    // is never used, which would have broken EVERY executable (clustering_hist,
    // rpt_v5_hist, export_training_data, ...) on the mu0/ttbar production where
    // those collections are dropped -- not just the rpt programs.
    //
    // Nothing lost: rpt_v5's paperIsPU labels a jet PU by "dR > 0.6 from all
    // truth HS jets with pT > 4 GeV", i.e. from truthHSJet* alone. The rpt_v2/3/4
    // references are comments only, and generate_rpt.cxx / scratch's
    // rpt_label_diag.cxx use the DIFFERENT per-jet
    // AntiKt4EMTopoJets_truthITPUJet_idx branches (neither is in the current
    // workflow). Re-add as std::unique_ptr members guarded on sample, the way
    // trackLeptonID is, if a future study actually needs them.

    TTreeReaderArray<std::vector<int>> topoJetTruthHSIdx;
    TTreeReaderArray<std::vector<int>> topoJetGhostTrackIdx;

    TTreeReaderArray<float> particleT;

    // Lepton–jet overlap removal (Z+jets only). Track_leptonID flags tracks that
    // are leptons; a non-zero value means "lepton" (works whether the branch
    // stores 0/1 or a signed PDG-style code). Bound lazily via unique_ptr only
    // when OVERLAP_REMOVAL is set, so samples whose ntuples lack the branch
    // (vbf/dijet/local) never touch it. removedJet is the per-event mask filled
    // by computeOverlapRemoval(); empty ⇒ nothing removed.
    // Track_leptonID is a `char` branch (holds the signed lepton PDG id: ±11/±13,
    // 0 for non-leptons — all within char's −128..127 range). It MUST be read
    // with a matching TTreeReaderArray<char>; a mismatched type (e.g. <int>)
    // silently yields size 0 (ROOT logs a CreateContentProxy error, which
    // clustering_hist's gErrorIgnoreLevel=kFatal hides), which would veto every
    // Z+jets event. Read individual values via leptonPdg() below so the sign is
    // correct regardless of whether `char` is signed on the build platform.
    std::unique_ptr<TTreeReaderArray<char>> trackLeptonID;
    std::vector<char> removedJet;

  BranchPointerWrapper(TTreeReader& r)
    : reader (r),
      weight (r, "weight"),
      trackD0 (r, "Track_d0"), trackZ0 (r, "Track_z0"), trackPt (r, "Track_pt"),
      trackEta (r, "Track_eta"), trackQP (r, "Track_qOverP"),
      trackTheta (r, "Track_theta"), trackPhi (r, "Track_phi"),
      trackVarD0 (r, "Track_var_d0"), trackVarZ0 (r, "Track_var_z0"),
      trackVarQp (r, "Track_var_qOverP"), trackVarTheta (r, "Track_var_theta"),
      trackVarPhi (r, "Track_var_phi0"),
      trackTime (r, "Track_time"), trackTimeRes (r, "Track_timeRes"),
      trackTimeValid (r, "Track_hasValidTime"), trackQuality (r, "Track_quality"),
      trackToTruthvtx (r, "Track_truthVtx_idx"),
      trackToParticle (r, "Track_truthPart_idx"),
      trackHgtdHits (r, "Track_nHGTDHits"),
      trackPrimHits (r, "Track_nHGTDPrimaryHits"),
      truthVtxZ (r, "TruthVtx_z"), truthVtxTime (r, "TruthVtx_time"),
      truthVtxIshs (r, "TruthVtx_isHS"), recoVtxZ (r, "RecoVtx_z"),
      recoVtxTime (r, "RecoVtx_time"), recoVtxTimeRes (r, "RecoVtx_timeRes"),
      recoVtxValid (r, "RecoVtx_hasValidTime"),
      topoJetPt (r, "AntiKt4EMTopoJets_pt"),
      topoJetEta (r, "AntiKt4EMTopoJets_eta"),
      topoJetPhi (r, "AntiKt4EMTopoJets_phi"),
      truthHSJetPt (r, "TruthHSJet_pt"), truthHSJetEta (r, "TruthHSJet_eta"),
      truthHSJetPhi (r, "TruthHSJet_phi"),
      topoJetTruthHSIdx (r, "AntiKt4EMTopoJets_truthHSJet_idx"),
      topoJetGhostTrackIdx (r, "AntiKt4EMTopoJets_ghostTrack_idx"),
      particleT (r, "TruthPart_prodVtx_time")
    {
      // Only bind the lepton branch when overlap removal is active (Z+jets),
      // so other samples' readers never register a branch their ntuples lack.
      if (OVERLAP_REMOVAL)
        trackLeptonID = std::make_unique<TTreeReaderArray<char>>(r, "Track_leptonID");
    }

    // -----------------------------------------------------------------------
    // leptonPdg
    //   Signed lepton PDG id for track t from the char Track_leptonID branch.
    //   The static_cast<signed char> reinterprets the stored 8-bit pattern as
    //   signed, so negative ids (e.g. -11) come back correctly even if `char`
    //   is unsigned on the build platform. Callers must ensure trackLeptonID is
    //   bound and t is in range.
    // -----------------------------------------------------------------------
    int leptonPdg(int t) const {
      return static_cast<int>(static_cast<signed char>((*trackLeptonID)[t]));
    }

    // -----------------------------------------------------------------------
    // isGoodLepton
    //   Track t is a selected lepton: Track_leptonID-flagged (non-zero signed
    //   PDG id) AND pt > LEPTON_MIN_PT. Single-sourced lepton definition shared
    //   by both the overlap removal and the Z-selection so the pT quality cut
    //   is applied consistently. Assumes t indexes a valid track (callers loop
    //   over the Track_leptonID array size).
    // -----------------------------------------------------------------------
    bool isGoodLepton(int t) const {
      return trackLeptonID && leptonPdg(t) != 0
             && this->trackPt[t] > LEPTON_MIN_PT;
    }

    // -----------------------------------------------------------------------
    // countGoodLeptons
    //   Number of tracks passing isGoodLepton in this event. Diagnostic-only
    //   helper (see LeptonSelDiag below) -- lets a low Z->ll survival rate be
    //   attributed to "too few qualifying leptons" vs. "found some but no
    //   OS-SF pair among them" instead of guessed at. Returns 0 when
    //   trackLeptonID isn't bound (non-Z+jets samples).
    // -----------------------------------------------------------------------
    int countGoodLeptons() const {
      if (!trackLeptonID) return 0;
      const int n = (int)trackLeptonID->GetSize();
      int count = 0;
      for (int t = 0; t < n; ++t)
        if (isGoodLepton(t)) ++count;
      return count;
    }

    // -----------------------------------------------------------------------
    // passLeptonSelection
    //   Z-boson event selection for Z+jets: require at least one opposite-sign
    //   same-flavour lepton pair among the selected leptons (isGoodLepton) —
    //   |pdg_i| == |pdg_j| and pdg_i·pdg_j < 0. Skips events with no
    //   reconstructed Z→ℓℓ. Returns true (no-op) when OVERLAP_REMOVAL is unset
    //   (vbf/dijet/local), leaving those samples unaffected. Unlike the
    //   SKIP_EVENT overlap veto, this is a clean physics selection (defines the
    //   Z→ℓℓ+jets fiducial region), not a pileup-biased cut.
    // -----------------------------------------------------------------------
    bool passLeptonSelection() const {
      if (!OVERLAP_REMOVAL || !trackLeptonID) return true;
      const int n = (int)trackLeptonID->GetSize();
      std::vector<int> ids;  // signed PDG ids of selected leptons
      for (int t = 0; t < n; ++t)
        if (isGoodLepton(t)) ids.push_back(leptonPdg(t));
      for (size_t a = 0; a < ids.size(); ++a)
        for (size_t b = a + 1; b < ids.size(); ++b)
          if (std::abs(ids[a]) == std::abs(ids[b]) && ids[a] * ids[b] < 0)
            return true;  // opposite-sign same-flavour pair found
      return false;
    }

    // -----------------------------------------------------------------------
    // computeOverlapRemoval
    //   Rebuilds removedJet for the current event: a reco jet is flagged
    //   removed if it lies within LEPTON_JET_DR of any selected lepton
    //   (isGoodLepton — flagged AND pt > LEPTON_MIN_PT; standard ATLAS
    //   lepton–jet overlap removal, all η). Must be called once per event
    //   before any jet loop. No-op (mask cleared to all-false) unless
    //   OVERLAP_REMOVAL is set, so on vbf/dijet/local isJetRemoved() is always
    //   false and behaviour is unchanged. Cheap: leptons per event are O(2).
    // -----------------------------------------------------------------------
    void computeOverlapRemoval() {
      if (!OVERLAP_REMOVAL || !trackLeptonID) { removedJet.clear(); return; }
      const int nJets = (int)this->topoJetPt.GetSize();
      removedJet.assign(nJets, 0);

      const int nLep = (int)trackLeptonID->GetSize();
      for (int t = 0; t < nLep; ++t) {
        if (!isGoodLepton(t)) continue;  // flagged lepton with pt > LEPTON_MIN_PT
        double lEta = this->trackEta[t];
        double lPhi = this->trackPhi[t];
        for (int j = 0; j < nJets; ++j) {
          if (removedJet[j]) continue;  // already removed by another lepton
          double deta = this->topoJetEta[j] - lEta;
          double dphi = TVector2::Phi_mpi_pi(this->topoJetPhi[j] - lPhi);
          if (std::hypot(deta, dphi) < LEPTON_JET_DR) removedJet[j] = 1;
        }
      }
    }

    // -----------------------------------------------------------------------
    // isJetRemoved
    //   True if reco jet j was dropped by lepton–jet overlap removal this event.
    //   Only meaningful in OverlapMode::REMOVE_JETS — in SKIP_EVENT mode the
    //   whole event is vetoed up front (see vetoLeptonOverlap), so no per-jet
    //   removal is applied and this returns false. Also always false when
    //   removedJet is empty (non-Z+jets samples, or before
    //   computeOverlapRemoval() has run), so every jet loop can guard on it
    //   unconditionally with zero effect on other samples.
    // -----------------------------------------------------------------------
    bool isJetRemoved(int j) const {
      if (OVERLAP_MODE != OverlapMode::REMOVE_JETS) return false;
      return j >= 0 && j < (int)removedJet.size() && removedJet[j] != 0;
    }

    // -----------------------------------------------------------------------
    // vetoLeptonOverlap
    //   True if this event should be skipped entirely because of a lepton–jet
    //   overlap. Only fires in OverlapMode::SKIP_EVENT (and only for Z+jets,
    //   via OVERLAP_REMOVAL); returns false in REMOVE_JETS mode, where the
    //   overlap is handled per-jet via isJetRemoved instead. Call after
    //   computeOverlapRemoval() has populated removedJet.
    // -----------------------------------------------------------------------
    bool vetoLeptonOverlap() const {
      if (!OVERLAP_REMOVAL || OVERLAP_MODE != OverlapMode::SKIP_EVENT) return false;
      for (char c : removedJet)
        if (c) return true;
      return false;
    }

    // -----------------------------------------------------------------------
    // passBasicCuts
    //   Returns false (skip event) if:
    //     • fewer than MIN_JETS reco jets are present (fast pre-filter before
    //       the full HS-matching check in passJetPtCut), or
    //     • the reconstructed HS vertex z is more than MAX_VTX_DZ from the
    //       truth HS vertex z (guards against wrong-vertex events).
    // -----------------------------------------------------------------------
    bool passBasicCuts() {
      if (this->topoJetPt.GetSize() < MIN_JETS) {
	if (DEBUG) std::cout << "Skipping low jet event\n";
        return false;
      }
    
      // check reco HS vertex is with 2mm of truth HS vertex
      if(std::abs(this->truthVtxZ[0] - this->recoVtxZ[0]) > MAX_VTX_DZ) {
	if(DEBUG) std::cout << "Skipping event due to incorrect HS vertex\n";
	return false;
      }
      
      return true;
    }

    // -----------------------------------------------------------------------
    // VbsPair / calcBestVbsPair
    //   The VBS candidate jet pair: among the given pT-passing jet indices, the
    //   pair in opposite η hemispheres (η_i · η_j < 0) with the largest
    //   invariant mass m_jj — the standard VBS-candidate-jet definition.
    //
    //   Returns m_jj alongside |Δη| (and the two reco jet indices) because the
    //   pair search is identical either way: this used to return |Δη| only and
    //   discard the m_jj it had just computed, which forced any m_jj study to
    //   redo the whole search. The indices let callers ask whether each leg is
    //   truth-HS-matched (isJetTruthHS) — i.e. whether the topology is genuine
    //   or manufactured by a pileup jet.
    //
    //   All fields keep their -1 sentinel when no opposite-hemisphere pair
    //   exists; passJetPtCut relies on dEta staying -1 there so a no-pair event
    //   fails even when the |Δη| requirement is disabled via --vbs-deta=0.
    // -----------------------------------------------------------------------
    struct VbsPair {
      double mjj  = -1.0;
      double dEta = -1.0;
      int    idxI = -1, idxJ = -1;
      bool valid() const { return idxI >= 0; }
    };

    VbsPair calcBestVbsPair(const std::vector<int>& passPtIdx) const {
      VbsPair best;
      for (size_t a = 0; a < passPtIdx.size(); ++a) {
        for (size_t b = a + 1; b < passPtIdx.size(); ++b) {
          int i = passPtIdx[a], j = passPtIdx[b];
          float etaI = this->topoJetEta[i], etaJ = this->topoJetEta[j];
          if (etaI * etaJ >= 0) continue;  // require opposite hemispheres

          TLorentzVector ji, jj;
          ji.SetPtEtaPhiM(this->topoJetPt[i], etaI, this->topoJetPhi[i], 0.0);
          jj.SetPtEtaPhiM(this->topoJetPt[j], etaJ, this->topoJetPhi[j], 0.0);
          double mjj = (ji + jj).M();

          if (mjj > best.mjj)
            best = {mjj, std::abs((double)etaI - (double)etaJ), i, j};
        }
      }
      return best;
    }

    // -----------------------------------------------------------------------
    // isJetTruthHS
    //   True if reco jet j is matched to at least one truth hard-scatter jet,
    //   via the ntuple's own AntiKt4EMTopoJets_truthHSJet_idx association (a
    //   non-empty index list means matched — the same idiom event_display.py
    //   uses). Preferred over the ΔR-cone paperIsHS lambdas duplicated across
    //   the rpt_v* diagnostics: it is the production's own matching, and the
    //   branch was already bound here but read by nothing.
    // -----------------------------------------------------------------------
    bool isJetTruthHS(int j) const {
      return j >= 0 && j < (int)this->topoJetTruthHSIdx.GetSize()
             && !this->topoJetTruthHSIdx[j].empty();
    }

    // -----------------------------------------------------------------------
    // isJetPaperHS / isJetPaperPU
    //   Jet truth labels per ATL-HGTD-PUB-2022-001 Sec. 3, hoisted out of the
    //   rpt_v* diagnostics (where they lived as event-loop-local lambdas,
    //   copy-pasted across rpt_v2/v3/v4/v5) so the clustering side can share
    //   them.
    //
    //   HS: dR < 0.3 from a truth HS jet with pT > 10 GeV.
    //   PU: dR > 0.6 from EVERY truth HS jet with pT > 4 GeV.
    //   Note PU is NOT !HS -- a jet in the 0.3-0.6 band, or near a 4-10 GeV
    //   truth jet, is neither, and callers are expected to skip those.
    //
    //   Deliberately kept separate from isJetTruthHS (the ntuple's own
    //   AntiKt4EMTopoJets_truthHSJet_idx association) rather than unified:
    //   classifyVbsRegion below uses THESE, because the R1/R2 ROCs in rpt_v5
    //   were measured with them, and switching definitions would silently
    //   redefine which events land in each region.
    // -----------------------------------------------------------------------
    bool isJetPaperHS(double jEta, double jPhi) const {
      for (int t = 0; t < (int)this->truthHSJetPt.GetSize(); ++t) {
        if (this->truthHSJetPt[t] < 10.0) continue;
        double deta = jEta - this->truthHSJetEta[t];
        double dphi = TVector2::Phi_mpi_pi(jPhi - this->truthHSJetPhi[t]);
        if (std::hypot(deta, dphi) < 0.3) return true;
      }
      return false;
    }

    bool isJetPaperPU(double jEta, double jPhi) const {
      for (int t = 0; t < (int)this->truthHSJetPt.GetSize(); ++t) {
        if (this->truthHSJetPt[t] < 4.0) continue;
        double deta = jEta - this->truthHSJetEta[t];
        double dphi = TVector2::Phi_mpi_pi(jPhi - this->truthHSJetPhi[t]);
        if (std::hypot(deta, dphi) < 0.6) return false;
      }
      return true;
    }

    // -----------------------------------------------------------------------
    // VbsRegion / classifyVbsRegion
    //   The three VBS topologies where forward timing bears on the event,
    //   classified off the VBS candidate pair (calcBestVbsPair):
    //
    //     R1  both legs forward (HGTD acceptance), one paper-HS + one paper-PU.
    //         Both jets are timeable, so timing must say WHICH is the
    //         hard-scatter one.
    //     R2  forward paper-PU leg + central paper-HS leg.
    //         Only the fake is timeable; the question is whether timing can
    //         reject it.
    //     R3  forward paper-HS leg + central paper-PU leg -- the mirror of R2.
    //         Only the GENUINE tag is timeable, so the question is whether a
    //         timing gate can accidentally empty it: there is no fake here for
    //         timing to reject, since the pileup leg is central and outside
    //         HGTD's reach regardless of what timing does.
    //
    //   Single source of truth shared by rpt_v5 (which fills per-jet RpT
    //   histograms per region) and the clustering analysis (whose VBF_R1 /
    //   VBF_R2 scores gate their denominator on the event's region). Both must
    //   agree on membership or the two "R1"s silently mean different samples.
    //
    //   The pair must also clear the runtime VBS topology knobs
    //   (VBS_JET_MJJ / VBS_JET_D_ETA) -- calcBestVbsPair only finds the
    //   max-m_jj pair, it does not cut on it.
    //
    //   outFwdHS / outFwdPU return the reco jet indices of the forward
    //   hard-scatter and forward pileup legs (-1 when the region has no such
    //   leg): R2 always has outFwdHS == -1 (its HS leg is central, not
    //   forward), R3 always has outFwdPU == -1 (its PU leg is central).
    // -----------------------------------------------------------------------
    VbsRegion classifyVbsRegion(double fwdEtaMin, double fwdEtaMax,
                                double centralEtaMax,
                                int* outFwdHS = nullptr,
                                int* outFwdPU = nullptr) const {
      if (outFwdHS) *outFwdHS = -1;
      if (outFwdPU) *outFwdPU = -1;

      std::vector<int> passPtIdx;
      int nPt = 0, nPtEta = 0;
      this->collectPtPassingJets(passPtIdx, nPt, nPtEta);
      VbsPair pair = this->calcBestVbsPair(passPtIdx);
      if (!pair.valid())              return VbsRegion::NONE;
      if (pair.mjj  < VBS_JET_MJJ)    return VbsRegion::NONE;
      if (pair.dEta < VBS_JET_D_ETA)  return VbsRegion::NONE;

      auto isFwd = [&](int j) {
        double e = std::abs((double)this->topoJetEta[j]);
        return e > fwdEtaMin && e < fwdEtaMax;
      };
      auto isCentral = [&](int j) {
        return std::abs((double)this->topoJetEta[j]) < centralEtaMax;
      };
      auto paperHS = [&](int j) {
        return this->isJetPaperHS(this->topoJetEta[j], this->topoJetPhi[j]);
      };
      auto paperPU = [&](int j) {
        return this->isJetPaperPU(this->topoJetEta[j], this->topoJetPhi[j]);
      };

      const int a = pair.idxI, b = pair.idxJ;

      // R1: both legs forward, exactly one HS and one PU.
      if (isFwd(a) && isFwd(b)) {
        if (paperHS(a) && paperPU(b)) {
          if (outFwdHS) *outFwdHS = a;
          if (outFwdPU) *outFwdPU = b;
          return VbsRegion::R1;
        }
        if (paperHS(b) && paperPU(a)) {
          if (outFwdHS) *outFwdHS = b;
          if (outFwdPU) *outFwdPU = a;
          return VbsRegion::R1;
        }
        return VbsRegion::NONE;
      }

      // R2: forward PU leg + central HS leg.
      if (isFwd(a) && paperPU(a) && isCentral(b) && paperHS(b)) {
        if (outFwdPU) *outFwdPU = a;
        return VbsRegion::R2;
      }
      if (isFwd(b) && paperPU(b) && isCentral(a) && paperHS(a)) {
        if (outFwdPU) *outFwdPU = b;
        return VbsRegion::R2;
      }

      // R3: forward HS leg + central PU leg. Written as R2's mirror on
      // purpose -- swap paperHS/paperPU and outFwdPU/outFwdHS -- so the
      // symmetry between "reject the fake" (R2) and "keep the genuine tag"
      // (R3) is visible in the code, not just asserted in the comment above.
      if (isFwd(a) && paperHS(a) && isCentral(b) && paperPU(b)) {
        if (outFwdHS) *outFwdHS = a;
        return VbsRegion::R3;
      }
      if (isFwd(b) && paperHS(b) && isCentral(a) && paperPU(a)) {
        if (outFwdHS) *outFwdHS = b;
        return VbsRegion::R3;
      }
      return VbsRegion::NONE;
    }

    // -----------------------------------------------------------------------
    // classifyEventRegion
    //   Event-level region membership -- see the EventRegion doc comment in
    //   clustering_constants.h for what each region means and why this exists
    //   alongside the pair-level classifyVbsRegion.
    //
    //   Uses the paper's dR-cone labels (isJetPaperHS / isJetPaperPU), the same
    //   definition the R_pT histogram fill uses, so region membership and the
    //   signal/background split inside a region cannot disagree.
    //
    //   Tried and rejected: isJetTruthHS (the ntuple's own truthHSJet_idx
    //   link). It is strictly binary, which looked attractive -- no "neither"
    //   jets -- but it tags 73.5% of forward jets above 30 GeV as hard-scatter
    //   against the paper definition's far more conservative labelling, so
    //   almost every event ended up with a forward HS jet and no forward
    //   pileup jet: 99.93% of events landed in a "HGTD can contribute" region,
    //   which is true under that definition and useless as a statement. NOT a
    //   soft-truth-jet artefact -- requiring the linked truth jet above 10 GeV
    //   moves it only to 72.6%, and above 30 GeV to 68.3%.
    //
    //   Consequence of the paper labels: they are NOT complements, so a jet
    //   can be neither HS nor PU (it contributes to no count below). Events
    //   still partition exactly across the four regions; individual jets need
    //   not.
    //
    //   Jets considered are exactly collectPtPassingJets' -- above MIN_JET_PT
    //   and not lepton-overlap-removed -- so region membership can never
    //   disagree with the jets the rest of the selection sees.
    //
    //   outCounts, when given, receives {nFwdHS, nFwdPU, nCenHS, nAnyHS} so a
    //   caller can report the sub-case breakdown without re-deriving it.
    // -----------------------------------------------------------------------
    EventRegion classifyEventRegion(double fwdEtaMin, double fwdEtaMax,
                                    double centralEtaMax,
                                    int* outCounts = nullptr) const {
      std::vector<int> passPtIdx;
      int nPt = 0, nPtEta = 0;
      this->collectPtPassingJets(passPtIdx, nPt, nPtEta);

      int nFwdHS = 0, nFwdPU = 0, nCenHS = 0, nAnyHS = 0;
      for (int j : passPtIdx) {
        const double e   = std::abs((double)this->topoJetEta[j]);
        const double phi = this->topoJetPhi[j];
        const double eta = this->topoJetEta[j];
        const bool   hs  = this->isJetPaperHS(eta, phi);
        const bool   pu  = this->isJetPaperPU(eta, phi);
        const bool  fwd  = (e > fwdEtaMin && e < fwdEtaMax);
        const bool  cen  = (e < centralEtaMax);
        if (hs) {
          ++nAnyHS;
          if (fwd) ++nFwdHS;
          else if (cen) ++nCenHS;
        } else if (pu && fwd) {
          ++nFwdPU;
        }
      }
      if (outCounts) {
        outCounts[0] = nFwdHS; outCounts[1] = nFwdPU;
        outCounts[2] = nCenHS; outCounts[3] = nAnyHS;
      }

      // Every reachable region requires a forward PU jet -- something for
      // timing to actually act against. A forward HS jet with NO forward
      // pileup is deliberately NOT reachable: timing would only confirm a tag
      // that had no competitor to be confused with, and on local VBF that
      // single case is 72.7% of events, which alone drove the "HGTD can
      // contribute" fraction to 98.9% and made the split meaningless.
      if (nFwdPU >= 1) {
        if (nFwdHS >= 1) return EventRegion::R1;        // both timeable
        if (nCenHS >= 1) return EventRegion::R2;        // only the fake timeable
        if (nAnyHS == 0) return EventRegion::CAN_HELP;  // fake, nothing genuine
        // else: the only hard-scatter jet sits beyond fwdEtaMax -- falls through.
      }
      return EventRegion::MAY_NOT;
    }

    // -----------------------------------------------------------------------
    // classifyR3Broad
    //   Three more cases -- on TOP of strict R3 -- where a timing gate still
    //   has a meaningful job to do, only ever called for a pair
    //   classifyVbsRegion has already returned NONE for (this function does
    //   NOT re-derive strict R1/R2/R3; callers add those separately from
    //   classifyVbsRegion's own return value, so a pair is never matched by
    //   both):
    //
    //     fwd HS + PU beyond fwdEtaMax   R3's story (protect the genuine tag)
    //                                    with the untimeable leg past the
    //                                    forward window instead of central.
    //     fwd PU + HS beyond fwdEtaMax   R2's story (reject the fake) with the
    //                                    untimeable leg past the window
    //                                    instead of central.
    //     both legs paper-PU, >=1 fwd    no genuine tag exists at all, so
    //                                    every forward (timeable) leg is a
    //                                    rejection candidate on its own --
    //                                    fills BOTH outPuLegA and outPuLegB
    //                                    when both legs happen to be forward.
    //
    //   At most one of {outHsLeg} or {outPuLegA[, outPuLegB]} is ever set per
    //   event: a pair is either a "confirm the genuine tag" case or a
    //   "reject the fake(s)" case, never both. Grouped under ONE name (R3)
    //   because that is how it is presented downstream -- one composition-plot
    //   band, one R_pT scenario set -- not because the two questions are the
    //   same one.
    // -----------------------------------------------------------------------
    bool classifyR3Broad(double fwdEtaMin, double fwdEtaMax, double centralEtaMax,
                         int* outHsLeg = nullptr,
                         int* outPuLegA = nullptr,
                         int* outPuLegB = nullptr) const {
      if (outHsLeg)  *outHsLeg  = -1;
      if (outPuLegA) *outPuLegA = -1;
      if (outPuLegB) *outPuLegB = -1;

      std::vector<int> passPtIdx;
      int nPt = 0, nPtEta = 0;
      this->collectPtPassingJets(passPtIdx, nPt, nPtEta);
      VbsPair pair = this->calcBestVbsPair(passPtIdx);
      if (!pair.valid())              return false;
      if (pair.mjj  < VBS_JET_MJJ)    return false;
      if (pair.dEta < VBS_JET_D_ETA)  return false;

      // Three-way zone, unlike classifyVbsRegion's two (isFwd/isCentral):
      // 0 = central, 1 = forward (HGTD acceptance), 2 = beyond fwdEtaMax. The
      // "beyond" zone is exactly the gap strict R2/R3 leave uncovered, since
      // both only ever pair a forward leg with a CENTRAL one.
      auto zone = [&](int j) {
        double e = std::abs((double)this->topoJetEta[j]);
        if (e < centralEtaMax)                    return 0;
        if (e > fwdEtaMin && e < fwdEtaMax)        return 1;
        return 2;
      };
      auto paperHS = [&](int j) {
        return this->isJetPaperHS(this->topoJetEta[j], this->topoJetPhi[j]);
      };
      auto paperPU = [&](int j) {
        return this->isJetPaperPU(this->topoJetEta[j], this->topoJetPhi[j]);
      };

      const int a = pair.idxI, b = pair.idxJ;
      const int zA = zone(a), zB = zone(b);
      const bool hsA = paperHS(a), puA = paperPU(a);
      const bool hsB = paperHS(b), puB = paperPU(b);

      // fwd HS + PU beyond: confirm the genuine tag.
      if (zA == 1 && hsA && zB == 2 && puB) { if (outHsLeg) *outHsLeg = a; return true; }
      if (zB == 1 && hsB && zA == 2 && puA) { if (outHsLeg) *outHsLeg = b; return true; }

      // fwd PU + HS beyond: reject the fake.
      if (zA == 1 && puA && zB == 2 && hsB) { if (outPuLegA) *outPuLegA = a; return true; }
      if (zB == 1 && puB && zA == 2 && hsA) { if (outPuLegA) *outPuLegA = b; return true; }

      // Both legs pileup: reject whichever leg(s) are actually timeable.
      if (puA && puB) {
        int* slot[2] = {outPuLegA, outPuLegB};
        int idx = 0;
        bool any = false;
        if (zA == 1) { if (slot[idx]) *slot[idx] = a; ++idx; any = true; }
        if (zB == 1) { if (idx < 2 && slot[idx]) *slot[idx] = b; ++idx; any = true; }
        return any;
      }

      return false;
    }

    // -----------------------------------------------------------------------
    // collectPtPassingJets
    //   Single jet loop behind passJetPtCut: fills passPtIdx with the indices
    //   of reco jets above MIN_JET_PT (skipping lepton-overlap-removed ones)
    //   and reports how many of those are also in the forward HGTD acceptance.
    //   Split out so a diagnostic can reach the same jet counts and the same
    //   VBS candidate pair *without* the |Δη| requirement bundled in — i.e.
    //   see the m_jj/|Δη| distribution of everything the rest of the selection
    //   admits, which is what a cut on either has to be chosen from.
    // -----------------------------------------------------------------------
    void collectPtPassingJets(std::vector<int>& passPtIdx,
                              int& nPassPt, int& nPassPtEta) const {
      passPtIdx.clear();
      nPassPt = 0; nPassPtEta = 0;
      for (int jetIdx = 0; jetIdx < (int)this->topoJetPt.GetSize(); ++jetIdx) {
        if (this->isJetRemoved(jetIdx)) continue;  // lepton-overlap removed (Z+jets)
        float eta = std::abs(this->topoJetEta[jetIdx]);
        float pt  = this->topoJetPt[jetIdx];
        if (DEBUG) std::cout << "reco jet pt, eta: " << pt << ", " << eta << '\n';
        bool passpt    = pt > MIN_JET_PT;
        bool passpteta = passpt && eta > MIN_ABS_ETA_JET && eta < MAX_ABS_ETA_JET;
        nPassPt    += passpt    ? 1 : 0;
        nPassPtEta += passpteta ? 1 : 0;
        if (passpt) passPtIdx.push_back(jetIdx);
      }
    }

    // -----------------------------------------------------------------------
    // passJetPtCut
    //   Requires at least MIN_PASSPT_JETS reco jets above MIN_JET_PT,
    //   at least MIN_PASSETA_JETS of which are in the forward HGTD acceptance,
    //   and that the VBS candidate pair — the opposite-hemisphere,
    //   pT-passing jet pair with the largest m_jj — has an η separation ≥
    //   VBS_JET_D_ETA, and m_jj >= VBS_JET_MJJ (0 by default -- a no-op, since
    //   any valid pair has mjj > 0; see the comment on VBS_JET_MJJ). No
    //   truth-HS matching required.
    // -----------------------------------------------------------------------
    bool passJetPtCut() {
      int passptcount = 0, passptetacount = 0;
      std::vector<int> passPtIdx;  // indices of pT-passing jets
      collectPtPassingJets(passPtIdx, passptcount, passptetacount);

      bool passesPt   = passptcount   >= MIN_PASSPT_JETS;
      bool passesEta  = passptetacount >= MIN_PASSETA_JETS;
      VbsPair pair    = calcBestVbsPair(passPtIdx);
      bool passesDEta = pair.dEta >= VBS_JET_D_ETA;
      bool passesMjj  = pair.mjj  >= VBS_JET_MJJ;
      return passesPt && passesEta && passesDEta && passesMjj;
    }

    // -----------------------------------------------------------------------
    // passForwardHsTracks  [testing only]
    //   Gates on the number of forward hard-scatter tracks being within the
    //   window [MIN_NHS_TRACK, MAX_NHS_TRACK].  Used during development to
    //   restrict studies to events with a controlled HS track multiplicity;
    //   not applied in the primary analysis.
    // -----------------------------------------------------------------------
    bool passForwardHsTracks(int& nForwardHSTrack) {
      return
	MAX_NHS_TRACK >= nForwardHSTrack and
	nForwardHSTrack >= MIN_NHS_TRACK;
    }

    // -----------------------------------------------------------------------
    // countForwardJets
    //   Counts reco jets that fall in the HGTD η acceptance
    //   (MIN_ABS_ETA_JET < |η| < MAX_ABS_ETA_JET) and exceed MIN_JET_PT.
    //   The result is written into nForwardJet and used as the x-axis
    //   variable in several efficiency/resolution plots.
    // -----------------------------------------------------------------------
    void countForwardJets(int& nForwardJet) {
    nForwardJet = 0;
    for(int jetIdx = 0; jetIdx < this->topoJetEta.GetSize(); ++jetIdx) {
      if (this->isJetRemoved(jetIdx)) continue;  // lepton-overlap removed (Z+jets)
      float
	jetEta = std::abs(this->topoJetEta[jetIdx]),
	jetPt = this->topoJetPt[jetIdx];
      if (jetEta > MIN_ABS_ETA_JET and jetEta < MAX_ABS_ETA_JET and jetPt > MIN_JET_PT)
	nForwardJet++;
    }
  }

    // -----------------------------------------------------------------------
    // countTruthHSJets
    //   Counts truth hard-scatter jets in the WHOLE event (no η restriction)
    //   above MIN_JET_PT.  Used as the binning x-axis variable in place of the
    //   forward reco-jet count.
    // -----------------------------------------------------------------------
    void countTruthHSJets(int& nTruthHSJet) {
    nTruthHSJet = 0;
    for (int jetIdx = 0; jetIdx < (int)this->truthHSJetPt.GetSize(); ++jetIdx) {
      if (this->truthHSJetPt[jetIdx] > MIN_JET_PT)
	nTruthHSJet++;
    }
  }

    // -----------------------------------------------------------------------
    // countForwardTracks
    //   Splits a pre-selected track list into forward HS and PU components.
    //   A track is counted if it passes the HGTD η window, pT bounds, track
    //   quality flag, and (when checkTimeValid is true) has a valid HGTD time.
    //   Writes totals into nFTrack, nFTrackHS, nFTrackPU; used as x-axis
    //   variables in efficiency/resolution plots.
    // -----------------------------------------------------------------------
    void countForwardTracks(
      int& nFTrack, int& nFTrackHS, int& nFTrackPU,
      const std::vector<int>& tracks, bool checkTimeValid
    ) {
      nFTrack = 0; nFTrackHS = 0; nFTrackPU = 0;
      for(auto trkIdx: tracks) {
	double eta = std::abs(this->trackEta[trkIdx]);
	double pt = this->trackPt[trkIdx];
	bool quality = this->trackQuality[trkIdx] == true;
	bool hasValidTime = checkTimeValid ? (this->trackTimeValid[trkIdx] == 1) : true;
	// already know these pass association
	if (eta > MIN_ABS_ETA_TRACK and eta < MAX_ABS_ETA_TRACK and
	    pt > MIN_TRACK_PT_COUNT and pt < MAX_TRACK_PT and
	    quality and hasValidTime) {
	  nFTrack++;
	  // Fetch truth vertex index only for tracks that pass the gate
	  int truthVtx = this->trackToTruthvtx[trkIdx];
	  if (truthVtx != -1 and this->truthVtxIshs[truthVtx])
	    nFTrackHS++;
	  else
	    nFTrackPU++;
	}
      }
    }

    // -----------------------------------------------------------------------
    // calcJetptDRScore
    //   Returns jet_pT * exp(−ΔR) where ΔR is the angular distance to the
    //   nearest reco jet.  Higher values indicate the track is close to a
    //   high-pT jet; used as a per-track weight in FILTJET-style scoring.
    // -----------------------------------------------------------------------
    double calcJetptDRScore(int trkIdx) {
      double
	trkEta = this->trackEta[trkIdx],
	trkPhi = this->trackPhi[trkIdx];

      double minDR = 1e6;
      int minIdx = -1;
      for (int jetIdx=0; jetIdx < this->topoJetEta.GetSize(); ++jetIdx) {
	double
	  jetEta = this->topoJetEta[jetIdx],
	  jetPhi = this->topoJetPhi[jetIdx];
	double
	  deta = jetEta-trkEta,
	  dphi = TVector2::Phi_mpi_pi(jetPhi - trkPhi);
	double thisDR = std::hypot(deta, dphi);
	if (thisDR < minDR) {
	  minDR = thisDR;
	  minIdx = jetIdx;
	}
      }
      double returnScore = this->topoJetPt[minIdx]*std::exp(-minDR);
      return returnScore;
    }

    // -----------------------------------------------------------------------
    // calcTrkptDRScore
    //   Returns track_pT * exp(−ΔR) where ΔR is the angular distance to the
    //   nearest reco jet.  Analogous to calcJetptDRScore but uses the track's
    //   own pT rather than the jet pT; used as the per-track contribution to
    //   the TRKPTZ cluster score.
    // -----------------------------------------------------------------------
    double calcTrkptDRScore(int trkIdx) {
      double
	trkEta = this->trackEta[trkIdx],
	trkPhi = this->trackPhi[trkIdx];


      double minDR = 1e6;
      for (int jetIdx = 0; this->topoJetEta.GetSize() > jetIdx; ++jetIdx) {
        double
	  jetEta = this->topoJetEta[jetIdx],
	  jetPhi = this->topoJetPhi[jetIdx];

	double
	  deta = jetEta-trkEta,
	  dphi = TVector2::Phi_mpi_pi(jetPhi - trkPhi);

	double thisDR = std::sqrt(deta*deta + dphi*dphi);
	if (thisDR < minDR) {
	  minDR = thisDR;
	}
      }
      double returnScore = this->trackPt[trkIdx]*std::exp(-minDR);
      return returnScore;
    }

  };

  // ---------------------------------------------------------------------------
  // EventCounts
  //   Holds all per-event counts and pre-folded histogram fill values
  //   derived from a 3σ track selection.  Computed once at the top of
  //   processEventData and passed to the histogram fill helpers below.
  // ---------------------------------------------------------------------------
  struct EventCounts {
    // Raw forward track/jet counts
    int nForwardJet      = 0;
    int nTruthHSJet      = 0;   // truth HS jets in the whole event (> MIN_JET_PT)
    int nForwardTrack    = 0;
    int nForwardTrackHS  = 0;
    int nForwardTrackPU  = 0;
    double puRatio       = 0.0;
    // Pre-folded x-values for efficiency / purity histograms
    int    effFillValFjet;
    int    effFillValHSJet;
    int    effFillValTrack;
    int    effFillValHSTrack;
    int    effFillValPUTrack;
    double effFillValPURatio;

    EventCounts(BranchPointerWrapper* branch,
                const std::vector<int>& tracks,
                bool checkValidTimes) {
      branch->countForwardJets(nForwardJet);
      branch->countTruthHSJets(nTruthHSJet);
      branch->countForwardTracks(nForwardTrack, nForwardTrackHS, nForwardTrackPU,
                                 tracks, checkValidTimes);
      puRatio          = (double)nForwardTrackPU / (double)nForwardTrack;
      effFillValFjet    = folded(nForwardJet,     (int)FOLD_FJET);
      effFillValHSJet   = folded(nTruthHSJet,     (int)FOLD_FJET);
      effFillValTrack   = folded(nForwardTrack,   (int)FOLD_TRACK);
      effFillValHSTrack = folded(nForwardTrackHS, (int)FOLD_HS_TRACK);
      effFillValPUTrack = folded(nForwardTrackPU, (int)FOLD_PU_TRACK);
      effFillValPURatio = folded(puRatio,         FOLD_PU_FRAC);
    }
  };

  // ---------------------------------------------------------------------------
  // Cluster
  //   Represents a single time(-z₀) cluster produced by doSimultaneousClustering
  //   or doConeClustering.  The weighted-mean position is stored in values[],
  //   with per-dimension uncertainties in sigmas[].  allTimes holds the raw
  //   track times for spread diagnostics.  trackIndices lists the original
  //   track indices from the branch arrays so that truth-matching and feature
  //   extraction can look up any per-track quantity.
  //
  //   Member methods:
  //     operator== / !=  — equality by (values[0], sigmas[0], nConstituents)
  //     calcPurity        — fraction of cluster pT belonging to truth HS vertex
  //     updateScores      — computes TRKPTZ, WAVES, JET_T_REFINED, WAVES_MISCL/MISAS, TEST_MISAS scores
  //     timeSpread        — std-dev of raw track times within the cluster
  //     zSpread           — std-dev of track z₀ values within the cluster
  //     passEfficiency    — true if cluster time is within 3·PASS_SIGMA of truth
  //     calcFeatures      — extracts and normalises the 10 ML input features
  // ---------------------------------------------------------------------------
  struct Cluster {
    std::vector<double> values;
    std::vector<double> sigmas;
    std::vector<double> allTimes;
    std::vector<int> trackIndices;
    std::unordered_map<int,double> scores;
    double purity = 0.0;
    bool wasMerged = false;
    int nConstituents=1;

    // -----------------------------------------------------------------------
    // operator== / operator!=
    //   Two clusters are considered identical if they share the same
    //   weighted-mean time (values[0]), time uncertainty (sigmas[0]), and
    //   constituent count.  This is sufficient for the duplicate-removal
    //   checks inside chooseCluster.
    // -----------------------------------------------------------------------
    bool operator==(const Cluster& other) {
      bool sameValues = values.at(0) == other.values.at(0);
      bool sameSigmas = sigmas.at(0) == other.sigmas.at(0);
      bool sameConsts = nConstituents == other.nConstituents;
      // this SHOULD be sufficient
      return sameConsts and sameSigmas and sameValues;
    }

    bool operator!=(const Cluster& other) {
      return !(*this == other);
    }

    // -----------------------------------------------------------------------
    // calcPurity
    //   Computes the pT-weighted fraction of cluster tracks that originate
    //   from the hard-scatter vertex (truth vertex index 0).  Stores the
    //   result in this->purity.  Called only when a purity-gated score is
    //   active (calcPurityFlag == true) to avoid unnecessary truth-matching work.
    // -----------------------------------------------------------------------
    void calcPurity(BranchPointerWrapper *branch) {
      double num = 0.0, denom = 0.0;
      for (auto trk: this->trackIndices) {
	if (branch->trackToTruthvtx[trk] == 0) num += branch->trackPt[trk];
	denom += branch->trackPt[trk];
      }
      this->purity = num/denom;
    }

    // -----------------------------------------------------------------------
    // updateScores
    //   Populates the score map for all derived scoring algorithms after
    //   clustering is complete.  Called once per cluster in
    //   clusterTracksInTime.
    //
    //   TRKPTZ   — TRKPT weighted by exp(−1.5 · |Δz|), where Δz is the
    //              distance from the cluster z₀ to the reco primary vertex.
    //              When z₀ is available as values[1] it is used directly;
    //              otherwise a precision-weighted average over track z₀ is
    //              computed on the fly.
    // -----------------------------------------------------------------------
    void updateScores(
      BranchPointerWrapper *branch
    ) {
      // Call calcFeatures once — it returns the normalised feature vector and
      // the raw deltaZ (cluster z − reco vertex z).  Reusing the returned
      // deltaZ avoids a second precision-weighted z-average pass over tracks.
      auto [features, rawDeltaZ, rawDeltaZResunits] = this->calcFeatures(branch);

      // TRKPTZ: score = TRKPT × exp(−1.5 · |Δz|)
      if (this->values.size() > 1) {
        // Clustering was done with z₀ as a second dimension: use that value.
        double dz = std::abs(this->values.at(1) - branch->recoVtxZ[0]);
        this->scores[Score::TRKPTZ.id] =
          this->scores.at(Score::TRKPT.id) * std::exp(-1.5 * dz);
      } else {
        // Common case (usez0=false): reuse deltaZ already computed by calcFeatures.
        this->scores[Score::TRKPTZ.id] =
          this->scores.at(Score::TRKPT.id) * std::exp(-1.5 * std::abs(rawDeltaZ));
      }

      // WAVES: WAVeS-style score — Σ_i pT_i × pT_jet(i) / max(ΔR_i, DR_FLOOR)
      // multiplied by exp(−1.5|Δz_cluster|), where Δz is the pT-weighted cluster z centroid
      // minus the reco vertex z.  The cluster-level z-term is more effective than per-track
      // z-pull weighting because it averages over track-by-track z noise.
      // Linear pT (not squared) so a couple of high-pT time-misassigned tracks can't
      // outvote a larger time-coherent cluster.
      // Falls back to TRKPTZ if no qualifying forward jets exist.
      {
        double wavesSum = 0.0;
        for (int idx : this->trackIndices) {
          double trkEta  = branch->trackEta[idx];
          double trkPhi  = branch->trackPhi[idx];
          double trkPt   = branch->trackPt[idx];

          double minDR     = 1e6;
          double nearJetPt = 0.0;
          for (int j = 0; j < (int)branch->topoJetEta.GetSize(); ++j) {
            if (branch->isJetRemoved(j)) continue;  // lepton-overlap removed (Z+jets)
            double jEta = branch->topoJetEta[j];
            double jPt  = branch->topoJetPt[j];
            if (jPt < MIN_JET_PT) continue;
            if (std::abs(jEta) < MIN_ABS_ETA_JET || std::abs(jEta) > MAX_ABS_ETA_JET) continue;
            double deta = jEta - trkEta;
            double dphi = TVector2::Phi_mpi_pi(branch->topoJetPhi[j] - trkPhi);
            double dr   = std::hypot(deta, dphi);
            if (dr < minDR) { minDR = dr; nearJetPt = jPt; }
          }
          if (nearJetPt <= 0.0) continue;
          wavesSum += trkPt * nearJetPt
                      / std::max(minDR, WAVES_DR_FLOOR);
        }
        if (wavesSum > 0.0)
          this->scores[Score::WAVES.id] =
              wavesSum * std::exp(-1.5 * std::abs(rawDeltaZ));
        else
          this->scores[Score::WAVES.id] = this->scores.at(Score::TRKPTZ.id);
      }

      // JET_T_REFINED: dedicated collection (jet-filtered tracks at 2σ iterative);
      // cluster selected by TRKPTZ via the aux-collection path in selectClusters().
      this->scores[Score::JET_T_REFINED.id] = this->scores.at(Score::TRKPTZ.id);

      // WAVeS oracle variants: selected by the WAVeS score; denominator gates
      // (cluster purity / HS timing purity) applied at fill time in event_processing.h.
      this->scores[Score::WAVES_MISCL.id] = this->scores.at(Score::WAVES.id);
      this->scores[Score::WAVES_MISAS.id] = this->scores.at(Score::WAVES.id);
      // VBS region rows: same WAVeS selection, region gate applied at fill time.
      this->scores[Score::VBF_R1.id]      = this->scores.at(Score::WAVES.id);
      this->scores[Score::VBF_R2.id]      = this->scores.at(Score::WAVES.id);

      // TEST_MISAS uses TRKPTZ as its selection score; the purity gate is applied
      // at efficiency-check time in event_processing.h (both pass and total fills).
      this->scores[Score::TEST_MISAS.id] = this->scores.at(Score::TRKPTZ.id);
    }
    
    // -----------------------------------------------------------------------
    // timeSpread
    //   Returns the population standard deviation of the raw track times
    //   stored in allTimes.  Used as a diagnostic variable in error-analysis
    //   plots to characterise how well the cluster's tracks agree in time.
    // -----------------------------------------------------------------------
    double timeSpread() {
      // calculate stdev of times
      auto avg = std::accumulate(allTimes.begin(),allTimes.end(),0.0)/this->allTimes.size();
      auto ssd = 0.0;
      for (auto t: allTimes) {
	ssd += (avg - t)*(avg - t);
      }

      return std::sqrt(ssd/this->allTimes.size());
    }

    // -----------------------------------------------------------------------
    // zSpread
    //   Returns the population standard deviation of track z₀ values for
    //   the cluster constituents.  Analogous to timeSpread but in the
    //   longitudinal direction; used in error-analysis diagnostics.
    // -----------------------------------------------------------------------
    double zSpread(BranchPointerWrapper *bpw) {
      // calculate stdev of z0 values
      auto avg = 0.0;
      for (auto trk: trackIndices)
	avg += bpw->trackZ0[trk];
      avg *= 1./this->trackIndices.size();
      auto ssd = 0.0;
      for (auto t: trackIndices) {
	ssd += (avg - bpw->trackZ0[t])*(avg - bpw->trackZ0[t]);
      }
      return std::sqrt(ssd/this->trackIndices.size());
    }

    // -----------------------------------------------------------------------
    // passEfficiency
    //   Returns true if the cluster's weighted-mean time is within
    //   PASS_SIGMA of the truth hard-scatter vertex time.
    //   This is the primary efficiency criterion used when filling the
    //   pass/total histograms in event_processing.h.
    // -----------------------------------------------------------------------
    bool passEfficiency(BranchPointerWrapper *branch) const {
      if (DEBUG) std::cout << "Choosing pass score\n";
      if (this->values.size() == 0)
	return false;

      double diff = std::abs(this->values.at(0)-branch->truthVtxTime[0]);
      if (diff > PASS_SIGMA)
	return false;

      return true;
    }

    // -----------------------------------------------------------------------
    // calculateTime
    //   Returns the time to use for diff/resolution plots and passEfficiency
    //   for the given score.  Implemented out-of-line in event_processing.h
    //   (after passTrackVertexAssociation is defined) to avoid circular includes.
    // -----------------------------------------------------------------------
    double calculateTime(
        Score score,
	BranchPointerWrapper* branch,
	double idealRes = -1.0
    ) const;

    // -----------------------------------------------------------------------
    // calculatePurity
    //   Returns the purity to use for resolution and efficiency plots for the
    //   given score.  For most scores this is the pre-computed this->purity
    //   (fraction of cluster ΣpT from HS tracks).  For WAVES, purity
    //   is re-evaluated using only the constituent tracks that survive the
    //   dR < 0.4 HS-jet filter, mirroring the time computation in calculateTime.
    //   Implemented out-of-line in event_processing.h alongside calculateTime.
    // -----------------------------------------------------------------------
    double calculatePurity(
        Score score,
        BranchPointerWrapper* branch
    ) const;

    // -----------------------------------------------------------------------
    // calcFeatures
    //   Extracts the 8 features used as DNN input and applies stored
    //   normalisation (means/stds from normalization_params.json) to return a
    //   unit-normalised feature vector.  Features (in training order):
    //     0  delta_z              — cluster z relative to reco primary vertex
    //     1  delta_z_resunits     — delta_z in units of cluster_z_sigma
    //     2  cluster_z_sigma      — precision-weighted z₀ uncertainty
    //     3  cluster_d0           — precision-weighted transverse impact parameter
    //     4  cluster_d0_sigma     — d0 uncertainty
    //     5  cluster_qOverP       — precision-weighted charge-over-momentum
    //     6  cluster_qOverP_sigma — q/p uncertainty
    //     7  cluster_sumpt        — scalar sum of constituent track pT
    //
    //   Returns {normalised feature vector, raw deltaZ}.
    //   The raw deltaZ is returned alongside so updateScores can compute
    //   TRKPTZ without repeating the z weighted-average loop.
    // -----------------------------------------------------------------------
    std::tuple<std::vector<float>, float, float> calcFeatures(BranchPointerWrapper *branch) {
      float znum  = 0.f, zden  = 0.f;
      float dnum  = 0.f, dden  = 0.f;
      float qpnum = 0.f, qpden = 0.f;
      float sumpt = 0.f;

      for (auto trk: trackIndices) {
        float trkVarZ = branch->trackVarZ0[trk];
        float trkVarD = branch->trackVarD0[trk];
        float trkVarQ = branch->trackVarQp[trk];

        znum  += branch->trackZ0[trk] / trkVarZ;  zden  += 1.f / trkVarZ;
        dnum  += branch->trackD0[trk] / trkVarD;  dden  += 1.f / trkVarD;
        qpnum += branch->trackQP[trk] / trkVarQ;  qpden += 1.f / trkVarQ;
        sumpt += branch->trackPt[trk];
      }

      float clusterZ           = znum / zden;
      float clusterZSigma      = 1.0f / std::sqrt(zden);
      float clusterD0          = dnum / dden;
      float clusterD0Sigma     = 1.0f / std::sqrt(dden);
      float clusterQOverP      = qpnum / qpden;
      float clusterQOverPSigma = 1.0f / std::sqrt(qpden);

      float deltaZ         = clusterZ - branch->recoVtxZ[0];
      float deltaZResunits = deltaZ / clusterZSigma;

      // Normalization parameters from share/models/normalization_params.json
      static const float MEANS[8] = {
        0.005349026450717473,  // delta_z mean
        0.006244652963375047,  // delta_z_resunits mean
        0.47186453018432917,  // cluster_z_sigma mean
        -0.00025930208416971584,  // cluster_d0 mean
        0.041210234810337115,  // cluster_d0_sigma mean
        2.497333134513254e-07,  // cluster_qOverP mean
        2.399697181725376e-06,  // cluster_qOverP_sigma mean
        22.605882118982347  // cluster_sumpt mean
      };
      static const float STDS[8] = {
        1.2256897153543274,  // delta_z std
        1.8559666108975377,  // delta_z_resunits std
        0.3927622369867232,  // cluster_z_sigma std
        0.15349002263169267,  // cluster_d0 std
        0.025843144355708243,  // cluster_d0_sigma std
        2.2377624188670666e-05,  // cluster_qOverP std
        1.5046303393371763e-06,  // cluster_qOverP_sigma std
        23.948572592936106  // cluster_sumpt std
      };

      std::vector<float> features = {
        (deltaZ             - MEANS[0]) / STDS[0],
        (deltaZResunits     - MEANS[1]) / STDS[1],
        (clusterZSigma      - MEANS[2]) / STDS[2],
        (clusterD0          - MEANS[3]) / STDS[3],
        (clusterD0Sigma     - MEANS[4]) / STDS[4],
        (clusterQOverP      - MEANS[5]) / STDS[5],
        (clusterQOverPSigma - MEANS[6]) / STDS[6],
        (sumpt              - MEANS[7]) / STDS[7],
      };

      return {features, deltaZ, deltaZResunits};
    }
  };
  
}
#endif // STRUCTS_H
