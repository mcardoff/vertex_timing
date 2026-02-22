#!/usr/bin/env python3
"""
Analyze the effect of pileup removal based on closest reconstructed vertex.
This script applies event selection cuts and analyzes removed tracks.
"""

import numpy as np
import uproot
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import glob
from pathlib import Path

# Analysis constants (from clustering_constants.h)
MIN_JETS = 1
MIN_PASSPT_JETS = 0
MIN_PASSETA_JETS = 1
MIN_JETPT = 30.0
MIN_ABS_ETA_JET = 2.38
MAX_ABS_ETA_JET = 4.00
MIN_ABS_ETA_TRACK = 2.38
MAX_ABS_ETA_TRACK = 4.00
MIN_TRACK_PT = 1.0
MAX_TRACK_PT = 30.0
MAX_VTX_DZ = 5.0
MAX_NSIGMA = 3.0

def pass_basic_cuts(event):
    """Check basic event-level cuts."""
    # Check minimum jets
    if len(event['TruthHSJet_pt']) < MIN_JETS or len(event['AntiKt4EMTopoJets_pt']) < MIN_JETS:
        return False

    # Check HS vertex reconstruction quality
    if abs(event['TruthVtx_z'][0] - event['RecoVtx_z'][0]) > MAX_VTX_DZ:
        return False

    return True

def pass_jet_pt_cut(event):
    """Check jet pT and eta cuts."""
    passpt_count = 0
    passpteta_count = 0

    for i in range(len(event['AntiKt4EMTopoJets_eta'])):
        eta = abs(event['AntiKt4EMTopoJets_eta'][i])
        pt = event['TruthHSJet_pt'][i] if i < len(event['TruthHSJet_pt']) else 0

        if pt > MIN_JETPT:
            passpt_count += 1
            if MIN_ABS_ETA_JET < eta < MAX_ABS_ETA_JET:
                passpteta_count += 1

    return passpt_count >= MIN_PASSPT_JETS and passpteta_count >= MIN_PASSETA_JETS

def get_forward_tracks(event):
    """Get indices of tracks passing quality, kinematic, and timing cuts."""
    forward_tracks = []

    n_tracks = len(event['Track_z0'])
    for idx in range(n_tracks):
        # Quality cuts
        if not event['Track_quality'][idx]:
            continue

        # Timing validity cut
        if not event['Track_hasValidTime'][idx]:
            continue

        # Kinematic cuts
        eta = abs(event['Track_eta'][idx])
        pt = event['Track_pt'][idx]

        if not (MIN_ABS_ETA_TRACK < eta < MAX_ABS_ETA_TRACK):
            continue
        if not (MIN_TRACK_PT < pt < MAX_TRACK_PT):
            continue

        forward_tracks.append(idx)

    return forward_tracks

def analyze_pileup_removal(event, forward_tracks):
    """Apply pileup removal and analyze removed tracks.

    Only considers tracks within MAX_NSIGMA of the reconstructed primary vertex.
    """
    pu_removed_tracks = []
    pu_kept_tracks = []
    pu_rejected_3sigma = []  # Tracks rejected by 3sigma cut

    z_primary = event['RecoVtx_z'][0]

    for idx_trk in forward_tracks:
        z0_trk = event['Track_z0'][idx_trk]
        var_z0_trk = event['Track_var_z0'][idx_trk]

        if var_z0_trk <= 0:  # Skip invalid variance
            continue

        # First, check if track is within MAX_NSIGMA of primary vertex
        nsigma_to_primary = abs(z_primary - z0_trk) / np.sqrt(var_z0_trk)

        if nsigma_to_primary > MAX_NSIGMA:
            pu_rejected_3sigma.append(idx_trk)
            continue  # Skip tracks beyond 3sigma of primary

        # Now find closest vertex among ALL reconstructed vertices
        closest_nsigma = np.inf
        closest_reco_vtx = -1

        for idx_vtx in range(len(event['RecoVtx_z'])):
            z_vtx = event['RecoVtx_z'][idx_vtx]
            nsigma = abs(z_vtx - z0_trk) / np.sqrt(var_z0_trk)

            if nsigma < closest_nsigma:
                closest_nsigma = nsigma
                closest_reco_vtx = idx_vtx

        # Tracks NOT associated to primary vertex (idx=0) are removed
        if closest_reco_vtx != 0:
            pu_removed_tracks.append(idx_trk)
        else:
            pu_kept_tracks.append(idx_trk)

    return pu_removed_tracks, pu_kept_tracks, pu_rejected_3sigma

def collect_statistics(input_files, max_events=None):
    """Collect statistics across all events."""

    results = {
        # Removed tracks (by closest vertex != primary)
        'removed_truthvtx': [],
        'removed_pt': [],
        'removed_z0': [],
        'removed_eta': [],
        'removed_is_hs': [],
        'removed_nsigma_to_primary': [],
        'removed_closest_vtx_idx': [],

        # Kept tracks (closest to primary within 3sigma)
        'kept_truthvtx': [],
        'kept_pt': [],
        'kept_z0': [],
        'kept_eta': [],
        'kept_is_hs': [],
        'kept_nsigma_to_primary': [],

        # Rejected by 3sigma cut (beyond 3sigma of primary)
        'rejected_3sigma_truthvtx': [],
        'rejected_3sigma_pt': [],
        'rejected_3sigma_z0': [],
        'rejected_3sigma_eta': [],
        'rejected_3sigma_is_hs': [],
        'rejected_3sigma_nsigma_to_primary': [],

        # Event-level stats
        'n_forward_tracks': [],
        'n_removed': [],
        'n_kept': [],
        'n_rejected_3sigma': [],
        'n_reco_vertices': [],
    }

    n_events_processed = 0
    n_events_passed = 0

    for file_path in input_files:
        print(f"Processing: {file_path}")

        try:
            with uproot.open(file_path) as f:
                tree = f['ntuple']

                # Get all branches we need
                branches = [
                    'Track_z0', 'Track_var_z0', 'Track_pt', 'Track_eta',
                    'Track_quality', 'Track_truthVtx_idx', 'Track_hasValidTime',
                    'TruthVtx_z', 'TruthVtx_isHS',
                    'RecoVtx_z',
                    'TruthHSJet_pt', 'TruthHSJet_eta',
                    'AntiKt4EMTopoJets_pt', 'AntiKt4EMTopoJets_eta',
                ]

                # Process in batches for efficiency
                for batch in tree.iterate(branches, library='np', step_size='100 MB'):
                    n_events_in_batch = len(batch['Track_z0'])

                    for i in range(n_events_in_batch):
                        if max_events and n_events_processed >= max_events:
                            break

                        # Extract event data
                        event = {key: batch[key][i] for key in branches}

                        n_events_processed += 1
                        if n_events_processed % 1000 == 0:
                            print(f"  Processed {n_events_processed} events, {n_events_passed} passed cuts")

                        # Apply event selection
                        if not pass_basic_cuts(event):
                            continue
                        if not pass_jet_pt_cut(event):
                            continue

                        n_events_passed += 1

                        # Get forward tracks
                        forward_tracks = get_forward_tracks(event)

                        # Apply pileup removal
                        removed_tracks, kept_tracks, rejected_3sigma = analyze_pileup_removal(event, forward_tracks)

                        # Store event-level stats
                        results['n_forward_tracks'].append(len(forward_tracks))
                        results['n_removed'].append(len(removed_tracks))
                        results['n_kept'].append(len(kept_tracks))
                        results['n_rejected_3sigma'].append(len(rejected_3sigma))
                        results['n_reco_vertices'].append(len(event['RecoVtx_z']))

                        # Store removed track info
                        for idx in removed_tracks:
                            truthvtx_idx = event['Track_truthVtx_idx'][idx]
                            results['removed_truthvtx'].append(truthvtx_idx)
                            results['removed_pt'].append(event['Track_pt'][idx])
                            results['removed_z0'].append(event['Track_z0'][idx])
                            results['removed_eta'].append(event['Track_eta'][idx])

                            # Check if from HS
                            is_hs = False
                            if truthvtx_idx >= 0 and truthvtx_idx < len(event['TruthVtx_isHS']):
                                is_hs = event['TruthVtx_isHS'][truthvtx_idx]
                            results['removed_is_hs'].append(is_hs)

                            # Calculate nsigma to primary vertex
                            z0_trk = event['Track_z0'][idx]
                            var_z0_trk = event['Track_var_z0'][idx]
                            z_primary = event['RecoVtx_z'][0]
                            nsigma = abs(z_primary - z0_trk) / np.sqrt(var_z0_trk) if var_z0_trk > 0 else 999
                            results['removed_nsigma_to_primary'].append(nsigma)

                            # Find which vertex it was associated to
                            closest_vtx = -1
                            closest_nsig = np.inf
                            for vidx in range(len(event['RecoVtx_z'])):
                                ns = abs(event['RecoVtx_z'][vidx] - z0_trk) / np.sqrt(var_z0_trk) if var_z0_trk > 0 else 999
                                if ns < closest_nsig:
                                    closest_nsig = ns
                                    closest_vtx = vidx
                            results['removed_closest_vtx_idx'].append(closest_vtx)

                        # Store kept track info
                        for idx in kept_tracks:
                            truthvtx_idx = event['Track_truthVtx_idx'][idx]
                            results['kept_truthvtx'].append(truthvtx_idx)
                            results['kept_pt'].append(event['Track_pt'][idx])
                            results['kept_z0'].append(event['Track_z0'][idx])
                            results['kept_eta'].append(event['Track_eta'][idx])

                            # Check if from HS
                            is_hs = False
                            if truthvtx_idx >= 0 and truthvtx_idx < len(event['TruthVtx_isHS']):
                                is_hs = event['TruthVtx_isHS'][truthvtx_idx]
                            results['kept_is_hs'].append(is_hs)

                            # Calculate nsigma to primary vertex
                            z0_trk = event['Track_z0'][idx]
                            var_z0_trk = event['Track_var_z0'][idx]
                            z_primary = event['RecoVtx_z'][0]
                            nsigma = abs(z_primary - z0_trk) / np.sqrt(var_z0_trk) if var_z0_trk > 0 else 999
                            results['kept_nsigma_to_primary'].append(nsigma)

                        # Store rejected_3sigma track info
                        for idx in rejected_3sigma:
                            truthvtx_idx = event['Track_truthVtx_idx'][idx]
                            results['rejected_3sigma_truthvtx'].append(truthvtx_idx)
                            results['rejected_3sigma_pt'].append(event['Track_pt'][idx])
                            results['rejected_3sigma_z0'].append(event['Track_z0'][idx])
                            results['rejected_3sigma_eta'].append(event['Track_eta'][idx])

                            # Check if from HS
                            is_hs = False
                            if truthvtx_idx >= 0 and truthvtx_idx < len(event['TruthVtx_isHS']):
                                is_hs = event['TruthVtx_isHS'][truthvtx_idx]
                            results['rejected_3sigma_is_hs'].append(is_hs)

                            # Calculate nsigma to primary vertex
                            z0_trk = event['Track_z0'][idx]
                            var_z0_trk = event['Track_var_z0'][idx]
                            z_primary = event['RecoVtx_z'][0]
                            nsigma = abs(z_primary - z0_trk) / np.sqrt(var_z0_trk) if var_z0_trk > 0 else 999
                            results['rejected_3sigma_nsigma_to_primary'].append(nsigma)

                    if max_events and n_events_processed >= max_events:
                        break

        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            continue

        if max_events and n_events_processed >= max_events:
            break

    print(f"\nTotal events processed: {n_events_processed}")
    print(f"Events passing cuts: {n_events_passed}")

    # Convert to numpy arrays
    for key in results:
        results[key] = np.array(results[key])

    return results

def plot_results(results, output_pdf):
    """Create diagnostic plots."""

    with PdfPages(output_pdf) as pdf:
        # Page 1: Truth vertex association
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Truth Vertex Association for Removed Tracks', fontsize=14, fontweight='bold')

        # Truth vertex distribution
        ax = axes[0, 0]
        ax.hist(results['removed_truthvtx'], bins=50, alpha=0.7, color='red', label='Removed')
        ax.hist(results['kept_truthvtx'], bins=50, alpha=0.5, color='blue', label='Kept')
        ax.set_xlabel('Truth Vertex Index')
        ax.set_ylabel('Tracks')
        ax.set_yscale('log')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # HS vs PU breakdown
        ax = axes[0, 1]
        removed_hs = np.sum(results['removed_is_hs'])
        removed_pu = len(results['removed_is_hs']) - removed_hs
        kept_hs = np.sum(results['kept_is_hs'])
        kept_pu = len(results['kept_is_hs']) - kept_hs
        rejected_3sig_hs = np.sum(results['rejected_3sigma_is_hs'])
        rejected_3sig_pu = len(results['rejected_3sigma_is_hs']) - rejected_3sig_hs

        x = np.arange(2)
        width = 0.25
        ax.bar(x - width, [removed_hs, removed_pu], width, label='Removed (wrong vtx)', color='red', alpha=0.7)
        ax.bar(x, [kept_hs, kept_pu], width, label='Kept', color='blue', alpha=0.7)
        ax.bar(x + width, [rejected_3sig_hs, rejected_3sig_pu], width, label='Rejected (>3σ)', color='gray', alpha=0.7)
        ax.set_ylabel('Number of Tracks')
        ax.set_xticks(x)
        ax.set_xticklabels(['Hard Scatter', 'Pileup'])
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_yscale('log')

        # Efficiency metrics
        ax = axes[1, 0]
        total_hs = removed_hs + kept_hs + rejected_3sig_hs
        total_pu = removed_pu + kept_pu + rejected_3sig_pu
        total_within_3sig_hs = removed_hs + kept_hs
        total_within_3sig_pu = removed_pu + kept_pu
        hs_efficiency = kept_hs / total_within_3sig_hs * 100 if total_within_3sig_hs > 0 else 0
        pu_rejection = removed_pu / total_within_3sig_pu * 100 if total_within_3sig_pu > 0 else 0

        metrics_text = f"""
        WITHIN 3σ OF PRIMARY VERTEX
        ============================
        Total Tracks: {len(results['removed_is_hs']) + len(results['kept_is_hs'])}

        Hard Scatter Tracks:
          Total: {total_within_3sig_hs}
          Kept: {kept_hs} ({hs_efficiency:.1f}%)
          Removed: {removed_hs} ({100-hs_efficiency:.1f}%)

        Pileup Tracks:
          Total: {total_within_3sig_pu}
          Kept: {kept_pu} ({100-pu_rejection:.1f}%)
          Removed: {removed_pu} ({pu_rejection:.1f}%)

        REJECTED BY 3σ CUT
        ==================
        HS: {rejected_3sig_hs} ({rejected_3sig_hs/total_hs*100:.1f}% of all HS)
        PU: {rejected_3sig_pu} ({rejected_3sig_pu/total_pu*100:.1f}% of all PU)
        """
        ax.text(0.1, 0.5, metrics_text, transform=ax.transAxes,
                fontsize=10, verticalalignment='center', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax.axis('off')

        # Closest vertex index for removed tracks
        ax = axes[1, 1]
        ax.hist(results['removed_closest_vtx_idx'], bins=50, alpha=0.7, color='red')
        ax.set_xlabel('Closest RecoVtx Index')
        ax.set_ylabel('Removed Tracks')
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3)
        ax.set_title('Which vertex removed tracks were associated to')

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

        # Page 2: pT distributions
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Track pT Distributions', fontsize=14, fontweight='bold')

        # All tracks pT
        ax = axes[0, 0]
        ax.hist(results['removed_pt'], bins=50, range=(MIN_TRACK_PT, MAX_TRACK_PT),
                alpha=0.7, color='red', label='Removed', density=True)
        ax.hist(results['kept_pt'], bins=50, range=(MIN_TRACK_PT, MAX_TRACK_PT),
                alpha=0.5, color='blue', label='Kept', density=True)
        ax.set_xlabel('Track pT [GeV]')
        ax.set_ylabel('Normalized Counts')
        ax.set_yscale('log')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_title('All Tracks')

        # HS tracks pT
        ax = axes[0, 1]
        removed_hs_pt = results['removed_pt'][results['removed_is_hs']]
        kept_hs_pt = results['kept_pt'][results['kept_is_hs']]
        ax.hist(removed_hs_pt, bins=50, range=(MIN_TRACK_PT, MAX_TRACK_PT),
                alpha=0.7, color='red', label=f'Removed HS (N={len(removed_hs_pt)})', density=True)
        ax.hist(kept_hs_pt, bins=50, range=(MIN_TRACK_PT, MAX_TRACK_PT),
                alpha=0.5, color='blue', label=f'Kept HS (N={len(kept_hs_pt)})', density=True)
        ax.set_xlabel('Track pT [GeV]')
        ax.set_ylabel('Normalized Counts')
        ax.set_yscale('log')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_title('Hard Scatter Tracks Only')

        # PU tracks pT
        ax = axes[1, 0]
        removed_pu_pt = results['removed_pt'][~results['removed_is_hs']]
        kept_pu_pt = results['kept_pt'][~results['kept_is_hs']]
        ax.hist(removed_pu_pt, bins=50, range=(MIN_TRACK_PT, MAX_TRACK_PT),
                alpha=0.7, color='red', label=f'Removed PU (N={len(removed_pu_pt)})', density=True)
        ax.hist(kept_pu_pt, bins=50, range=(MIN_TRACK_PT, MAX_TRACK_PT),
                alpha=0.5, color='blue', label=f'Kept PU (N={len(kept_pu_pt)})', density=True)
        ax.set_xlabel('Track pT [GeV]')
        ax.set_ylabel('Normalized Counts')
        ax.set_yscale('log')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_title('Pileup Tracks Only')

        # Mean pT comparison
        ax = axes[1, 1]
        mean_values = [
            np.mean(removed_hs_pt) if len(removed_hs_pt) > 0 else 0,
            np.mean(kept_hs_pt) if len(kept_hs_pt) > 0 else 0,
            np.mean(removed_pu_pt) if len(removed_pu_pt) > 0 else 0,
            np.mean(kept_pu_pt) if len(kept_pu_pt) > 0 else 0,
        ]
        x = np.arange(4)
        colors = ['red', 'blue', 'red', 'blue']
        ax.bar(x, mean_values, color=colors, alpha=0.7)
        ax.set_ylabel('Mean pT [GeV]')
        ax.set_xticks(x)
        ax.set_xticklabels(['Removed\nHS', 'Kept\nHS', 'Removed\nPU', 'Kept\nPU'])
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_title('Mean Track pT')

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

        # Page 3: z0 distributions
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Track z0 Distributions', fontsize=14, fontweight='bold')

        # All tracks z0
        ax = axes[0, 0]
        ax.hist(results['removed_z0'], bins=100, range=(-200, 200),
                alpha=0.7, color='red', label='Removed', density=True)
        ax.hist(results['kept_z0'], bins=100, range=(-200, 200),
                alpha=0.5, color='blue', label='Kept', density=True)
        ax.set_xlabel('Track z0 [mm]')
        ax.set_ylabel('Normalized Counts')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_title('All Tracks')

        # HS tracks z0
        ax = axes[0, 1]
        removed_hs_z0 = results['removed_z0'][results['removed_is_hs']]
        kept_hs_z0 = results['kept_z0'][results['kept_is_hs']]
        ax.hist(removed_hs_z0, bins=100, range=(-200, 200),
                alpha=0.7, color='red', label=f'Removed HS (N={len(removed_hs_z0)})', density=True)
        ax.hist(kept_hs_z0, bins=100, range=(-200, 200),
                alpha=0.5, color='blue', label=f'Kept HS (N={len(kept_hs_z0)})', density=True)
        ax.set_xlabel('Track z0 [mm]')
        ax.set_ylabel('Normalized Counts')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_title('Hard Scatter Tracks Only')

        # PU tracks z0
        ax = axes[1, 0]
        removed_pu_z0 = results['removed_z0'][~results['removed_is_hs']]
        kept_pu_z0 = results['kept_z0'][~results['kept_is_hs']]
        ax.hist(removed_pu_z0, bins=100, range=(-200, 200),
                alpha=0.7, color='red', label=f'Removed PU (N={len(removed_pu_z0)})', density=True)
        ax.hist(kept_pu_z0, bins=100, range=(-200, 200),
                alpha=0.5, color='blue', label=f'Kept PU (N={len(kept_pu_z0)})', density=True)
        ax.set_xlabel('Track z0 [mm]')
        ax.set_ylabel('Normalized Counts')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_title('Pileup Tracks Only')

        # z0 RMS comparison
        ax = axes[1, 1]
        rms_values = [
            np.std(removed_hs_z0) if len(removed_hs_z0) > 0 else 0,
            np.std(kept_hs_z0) if len(kept_hs_z0) > 0 else 0,
            np.std(removed_pu_z0) if len(removed_pu_z0) > 0 else 0,
            np.std(kept_pu_z0) if len(kept_pu_z0) > 0 else 0,
        ]
        x = np.arange(4)
        colors = ['red', 'blue', 'red', 'blue']
        ax.bar(x, rms_values, color=colors, alpha=0.7)
        ax.set_ylabel('z0 RMS [mm]')
        ax.set_xticks(x)
        ax.set_xticklabels(['Removed\nHS', 'Kept\nHS', 'Removed\nPU', 'Kept\nPU'])
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_title('Track z0 Spread')

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

        # Page 4: nsigma to primary vertex
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Track Distance to Primary Vertex (nσ)', fontsize=14, fontweight='bold')

        # All removed tracks nsigma
        ax = axes[0, 0]
        ax.hist(results['removed_nsigma_to_primary'], bins=60, range=(0, 3),
                alpha=0.7, color='red', label='Removed')
        ax.hist(results['kept_nsigma_to_primary'], bins=60, range=(0, 3),
                alpha=0.5, color='blue', label='Kept')
        ax.axvline(MAX_NSIGMA, color='black', linestyle='--', linewidth=2, label=f'Cut threshold ({MAX_NSIGMA}σ)')
        ax.set_xlabel('nσ to Primary Vertex')
        ax.set_ylabel('Tracks')
        ax.set_yscale('log')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_title('All Tracks (Within 3σ)')
        ax.set_xlim(0, 3)

        # HS tracks nsigma
        ax = axes[0, 1]
        removed_hs_nsig = results['removed_nsigma_to_primary'][results['removed_is_hs']]
        kept_hs_nsig = results['kept_nsigma_to_primary'][results['kept_is_hs']]
        ax.hist(removed_hs_nsig, bins=60, range=(0, 3),
                alpha=0.7, color='red', label=f'Removed HS (N={len(removed_hs_nsig)})')
        ax.hist(kept_hs_nsig, bins=60, range=(0, 3),
                alpha=0.5, color='blue', label=f'Kept HS (N={len(kept_hs_nsig)})')
        ax.axvline(MAX_NSIGMA, color='black', linestyle='--', linewidth=2, label=f'Cut threshold ({MAX_NSIGMA}σ)')
        ax.set_xlabel('nσ to Primary Vertex')
        ax.set_ylabel('Tracks')
        ax.set_yscale('log')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_title('Hard Scatter Tracks (Within 3σ)')
        ax.set_xlim(0, 3)

        # PU tracks nsigma
        ax = axes[1, 0]
        removed_pu_nsig = results['removed_nsigma_to_primary'][~results['removed_is_hs']]
        kept_pu_nsig = results['kept_nsigma_to_primary'][~results['kept_is_hs']]
        ax.hist(removed_pu_nsig, bins=60, range=(0, 3),
                alpha=0.7, color='red', label=f'Removed PU (N={len(removed_pu_nsig)})')
        ax.hist(kept_pu_nsig, bins=60, range=(0, 3),
                alpha=0.5, color='blue', label=f'Kept PU (N={len(kept_pu_nsig)})')
        ax.axvline(MAX_NSIGMA, color='black', linestyle='--', linewidth=2, label=f'Cut threshold ({MAX_NSIGMA}σ)')
        ax.set_xlabel('nσ to Primary Vertex')
        ax.set_ylabel('Tracks')
        ax.set_yscale('log')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_title('Pileup Tracks (Within 3σ)')
        ax.set_xlim(0, 3)

        # Statistics
        ax = axes[1, 1]
        stats_text = f"""
        nσ to Primary Vertex Statistics
        (All tracks within 3σ by construction)
        =======================================

        Hard Scatter Tracks:
          Removed: mean={np.mean(removed_hs_nsig):.2f}, median={np.median(removed_hs_nsig):.2f}
          Kept:    mean={np.mean(kept_hs_nsig):.2f}, median={np.median(kept_hs_nsig):.2f}

        Pileup Tracks:
          Removed: mean={np.mean(removed_pu_nsig):.2f}, median={np.median(removed_pu_nsig):.2f}
          Kept:    mean={np.mean(kept_pu_nsig):.2f}, median={np.median(kept_pu_nsig):.2f}

        Fraction of removed HS tracks with nσ < 1:
          {np.sum(removed_hs_nsig < 1) / len(removed_hs_nsig) * 100:.1f}%

        Fraction of kept HS tracks with nσ < 1:
          {np.sum(kept_hs_nsig < 1) / len(kept_hs_nsig) * 100:.1f}%
        """
        ax.text(0.1, 0.5, stats_text, transform=ax.transAxes,
                fontsize=9, verticalalignment='center', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax.axis('off')

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

        # Page 5: Event-level statistics
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Event-Level Track Statistics', fontsize=14, fontweight='bold')

        # Forward tracks per event
        ax = axes[0, 0]
        ax.hist(results['n_forward_tracks'], bins=50, alpha=0.7, color='green')
        ax.set_xlabel('Number of Forward Tracks')
        ax.set_ylabel('Events')
        ax.grid(True, alpha=0.3)
        ax.set_title(f'Mean: {np.mean(results["n_forward_tracks"]):.1f}')

        # Removed vs kept
        ax = axes[0, 1]
        ax.scatter(results['n_removed'], results['n_kept'], alpha=0.1, s=1)
        ax.set_xlabel('Tracks Removed')
        ax.set_ylabel('Tracks Kept')
        ax.grid(True, alpha=0.3)
        ax.set_title('Per-Event Track Removal')

        # Fraction removed
        ax = axes[1, 0]
        frac_removed = results['n_removed'] / (results['n_removed'] + results['n_kept'])
        frac_removed = frac_removed[~np.isnan(frac_removed)]
        ax.hist(frac_removed, bins=50, range=(0, 1), alpha=0.7, color='orange')
        ax.set_xlabel('Fraction of Tracks Removed')
        ax.set_ylabel('Events')
        ax.grid(True, alpha=0.3)
        ax.set_title(f'Mean: {np.mean(frac_removed):.3f}')

        # Number of reconstructed vertices
        ax = axes[1, 1]
        ax.hist(results['n_reco_vertices'], bins=50, alpha=0.7, color='purple')
        ax.set_xlabel('Number of Reconstructed Vertices')
        ax.set_ylabel('Events')
        ax.grid(True, alpha=0.3)
        ax.set_title(f'Mean: {np.mean(results["n_reco_vertices"]):.1f}')

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

        # Page 6: Eta distributions
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Track η Distributions', fontsize=14, fontweight='bold')

        # All tracks eta
        ax = axes[0, 0]
        ax.hist(results['removed_eta'], bins=50, range=(MIN_ABS_ETA_TRACK, MAX_ABS_ETA_TRACK),
                alpha=0.7, color='red', label='Removed', density=True)
        ax.hist(results['kept_eta'], bins=50, range=(MIN_ABS_ETA_TRACK, MAX_ABS_ETA_TRACK),
                alpha=0.5, color='blue', label='Kept', density=True)
        ax.set_xlabel('Track |η|')
        ax.set_ylabel('Normalized Counts')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_title('All Tracks')

        # HS tracks eta
        ax = axes[0, 1]
        removed_hs_eta = results['removed_eta'][results['removed_is_hs']]
        kept_hs_eta = results['kept_eta'][results['kept_is_hs']]
        ax.hist(removed_hs_eta, bins=50, range=(MIN_ABS_ETA_TRACK, MAX_ABS_ETA_TRACK),
                alpha=0.7, color='red', label=f'Removed HS (N={len(removed_hs_eta)})', density=True)
        ax.hist(kept_hs_eta, bins=50, range=(MIN_ABS_ETA_TRACK, MAX_ABS_ETA_TRACK),
                alpha=0.5, color='blue', label=f'Kept HS (N={len(kept_hs_eta)})', density=True)
        ax.set_xlabel('Track |η|')
        ax.set_ylabel('Normalized Counts')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_title('Hard Scatter Tracks')

        # PU tracks eta
        ax = axes[1, 0]
        removed_pu_eta = results['removed_eta'][~results['removed_is_hs']]
        kept_pu_eta = results['kept_eta'][~results['kept_is_hs']]
        ax.hist(removed_pu_eta, bins=50, range=(MIN_ABS_ETA_TRACK, MAX_ABS_ETA_TRACK),
                alpha=0.7, color='red', label=f'Removed PU (N={len(removed_pu_eta)})', density=True)
        ax.hist(kept_pu_eta, bins=50, range=(MIN_ABS_ETA_TRACK, MAX_ABS_ETA_TRACK),
                alpha=0.5, color='blue', label=f'Kept PU (N={len(kept_pu_eta)})', density=True)
        ax.set_xlabel('Track |η|')
        ax.set_ylabel('Normalized Counts')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_title('Pileup Tracks')

        # Summary table
        ax = axes[1, 1]
        n_total_tracks = len(results['removed_is_hs']) + len(results['kept_is_hs']) + len(results['rejected_3sigma_is_hs'])
        summary_text = f"""
        SUMMARY STATISTICS
        ==================

        Total Events: {len(results['n_forward_tracks'])}

        All Forward Tracks: {n_total_tracks}
          - Within 3σ: {len(results['removed_is_hs']) + len(results['kept_is_hs'])}
            • Removed: {len(results['removed_is_hs'])}
            • Kept: {len(results['kept_is_hs'])}
          - Rejected (>3σ): {len(results['rejected_3sigma_is_hs'])}

        HS Efficiency (within 3σ): {hs_efficiency:.1f}%
        PU Rejection (within 3σ): {pu_rejection:.1f}%

        Avg tracks per event:
          - Forward: {np.mean(results['n_forward_tracks']):.1f}
          - Removed: {np.mean(results['n_removed']):.1f}
          - Kept: {np.mean(results['n_kept']):.1f}
          - Rejected (>3σ): {np.mean(results['n_rejected_3sigma']):.1f}

        Avg RecoVtx per event: {np.mean(results['n_reco_vertices']):.1f}
        """
        ax.text(0.1, 0.5, summary_text, transform=ax.transAxes,
                fontsize=10, verticalalignment='center', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))
        ax.axis('off')

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

        # Page 7: ROC Curves
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('ROC Curves: HS Efficiency vs PU Rejection', fontsize=14, fontweight='bold')

        # For ROC curves, we evaluate ONLY the closest-vertex decision
        # (the nσ pre-cut would be a separate analysis)
        # So we look at tracks that passed the 3σ cut and see how well closest-vertex works

        # All tracks that passed 3σ cut
        all_nsigma = np.concatenate([results['removed_nsigma_to_primary'], results['kept_nsigma_to_primary']])
        all_is_hs = np.concatenate([results['removed_is_hs'], results['kept_is_hs']])
        all_pt = np.concatenate([results['removed_pt'], results['kept_pt']])
        all_z0 = np.concatenate([results['removed_z0'], results['kept_z0']])

        # Track status: True = kept (associated to primary), False = removed (associated to other vtx)
        all_kept = np.concatenate([
            np.zeros(len(results['removed_is_hs']), dtype=bool),  # removed tracks
            np.ones(len(results['kept_is_hs']), dtype=bool)       # kept tracks
        ])

        # Calculate |delta_z| to primary for all tracks
        all_delta_z = np.abs(all_z0)

        # ROC 1: nσ to primary vertex as discriminant
        # Hypothesis: tracks with smaller nσ to primary are more likely HS
        ax = axes[0, 0]
        thresholds_nsigma = np.linspace(0, 3, 100)
        hs_eff_nsigma = []
        pu_rej_nsigma = []

        for thresh in thresholds_nsigma:
            # "Keep" tracks with nsigma < thresh (closer to primary)
            decision_keep = all_nsigma < thresh

            # HS efficiency: fraction of HS tracks we keep
            hs_kept = np.sum(decision_keep & all_is_hs)
            hs_total = np.sum(all_is_hs)
            hs_eff = hs_kept / hs_total * 100 if hs_total > 0 else 0

            # PU rejection: fraction of PU tracks we remove
            pu_removed = np.sum(~decision_keep & ~all_is_hs)
            pu_total = np.sum(~all_is_hs)
            pu_rej = pu_removed / pu_total * 100 if pu_total > 0 else 0

            hs_eff_nsigma.append(hs_eff)
            pu_rej_nsigma.append(pu_rej)

        ax.plot(hs_eff_nsigma, pu_rej_nsigma, 'b-', linewidth=2)
        # Note: current algorithm uses closest-vertex, not nσ threshold
        # So current point may not lie on this curve
        ax.scatter([hs_efficiency], [pu_rejection], color='red', s=100, zorder=5,
                   label=f'Current algo: ε={hs_efficiency:.1f}%, rej={pu_rejection:.1f}%', marker='*')
        ax.set_xlabel('HS Efficiency (%)')
        ax.set_ylabel('PU Rejection (%)')
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
        ax.set_title('ROC: nσ to Primary (simple threshold)')

        # ROC 2: |Δz| to primary vertex (scan 0-5mm)
        ax = axes[0, 1]
        thresholds_dz = np.linspace(0, 5, 100)
        hs_eff_dz = []
        pu_rej_dz = []

        for thresh in thresholds_dz:
            # Keep tracks with |delta_z| < thresh
            kept_mask = all_delta_z < thresh

            hs_kept = np.sum(kept_mask & all_is_hs)
            hs_total = np.sum(all_is_hs)
            hs_eff = hs_kept / hs_total * 100 if hs_total > 0 else 0

            pu_removed = np.sum(~kept_mask & ~all_is_hs)
            pu_total = np.sum(~all_is_hs)
            pu_rej = pu_removed / pu_total * 100 if pu_total > 0 else 0

            hs_eff_dz.append(hs_eff)
            pu_rej_dz.append(pu_rej)

        ax.plot(hs_eff_dz, pu_rej_dz, 'b-', linewidth=2)
        ax.set_xlabel('HS Efficiency (%)')
        ax.set_ylabel('PU Rejection (%)')
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        ax.grid(True, alpha=0.3)
        ax.set_title('ROC: |Δz| to Primary Vertex')

        # ROC 3: Track pT threshold (scan 1-30 GeV, INVERTED)
        ax = axes[1, 0]
        thresholds_pt = np.linspace(1, 30, 100)
        hs_eff_pt = []
        pu_rej_pt = []

        for thresh in thresholds_pt:
            # Keep tracks with pT > thresh (higher pT tracks)
            kept_mask = all_pt > thresh

            hs_kept = np.sum(kept_mask & all_is_hs)
            hs_total = np.sum(all_is_hs)
            hs_eff = hs_kept / hs_total * 100 if hs_total > 0 else 0

            pu_removed = np.sum(~kept_mask & ~all_is_hs)
            pu_total = np.sum(~all_is_hs)
            pu_rej = pu_removed / pu_total * 100 if pu_total > 0 else 0

            hs_eff_pt.append(hs_eff)
            pu_rej_pt.append(pu_rej)

        ax.plot(hs_eff_pt, pu_rej_pt, 'b-', linewidth=2)
        ax.set_xlabel('HS Efficiency (%)')
        ax.set_ylabel('PU Rejection (%)')
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        ax.grid(True, alpha=0.3)
        ax.set_title('ROC: Track pT > threshold')

        # ROC 4: Comparison plot with explanatory text
        ax = axes[1, 1]
        ax.plot(hs_eff_nsigma, pu_rej_nsigma, 'b-', linewidth=2, label='nσ threshold', alpha=0.7)
        ax.plot(hs_eff_dz, pu_rej_dz, 'g-', linewidth=2, label='|Δz| threshold', alpha=0.7)
        ax.plot(hs_eff_pt, pu_rej_pt, 'orange', linewidth=2, label='pT threshold', alpha=0.7)
        ax.scatter([hs_efficiency], [pu_rejection], color='red', s=150, zorder=5,
                   label=f'Current algo (closest vtx)', marker='*')
        # Add diagonal for reference
        ax.plot([0, 100], [0, 100], 'k--', alpha=0.3, linewidth=1, label='Random')

        # Add text box explaining the mismatch
        textstr = 'Note: Current algorithm uses\nclosest-vertex association,\nnot simple thresholds.\nRed star shows actual performance.'
        ax.text(0.05, 0.35, textstr, transform=ax.transAxes,
                fontsize=8, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        ax.set_xlabel('HS Efficiency (%)')
        ax.set_ylabel('PU Rejection (%)')
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7, loc='lower right')
        ax.set_title('ROC Comparison: Alternative Discriminants')

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

def main():
    import argparse

    parser = argparse.ArgumentParser(description='Analyze pileup removal process')
    parser.add_argument('--input', '-i', required=True, help='Input ROOT file(s) (supports wildcards)')
    parser.add_argument('--output', '-o', default='pileup_removal_analysis.pdf',
                        help='Output PDF file')
    parser.add_argument('--max-events', '-n', type=int, default=None,
                        help='Maximum number of events to process')

    args = parser.parse_args()

    # Find input files
    input_files = glob.glob(args.input)
    if not input_files:
        print(f"No files found matching: {args.input}")
        return

    print(f"Found {len(input_files)} input file(s)")

    # Collect statistics
    results = collect_statistics(input_files, max_events=args.max_events)

    # Create plots
    print(f"\nGenerating plots: {args.output}")
    plot_results(results, args.output)

    print("\nDone!")

if __name__ == '__main__':
    main()
