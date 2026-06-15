"""
EVENT DISPLAY CODE
CREATED BY: wasikul.islam@cern.ch
MODIFIED FOR USE BY: mcardiff@brandeis.edu
"""

import argparse
import random
import subprocess
import sys
from itertools import chain

import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

# --- Argument Parsing and Data Loading ---
parser = argparse.ArgumentParser(description='Process event and vertex data.')
parser.add_argument('--event_num', type=int, required=True)
parser.add_argument('--file_num', type=str, required=True)
parser.add_argument('--extra_time', type=float, required=False)
parser.add_argument('--jet_idx', type=int, required=False, default=None,
                    help='Index of jet to highlight (orange) in R-Z and eta-phi displays')
parser.add_argument('--jet_label', type=str, required=False, default=None,
                    help='Truth identity of the target jet: HS or PU')
parser.add_argument('--rpt_hgtd', type=float, required=False, default=None,
                    help='HGTD RpT value for the target jet')
parser.add_argument('--rpt_mine', type=float, required=False, default=None,
                    help='My algo RpT value for the target jet')
parser.add_argument('--output_dir', type=str, required=False,
                    default='event_displays',
                    help='Directory to save the output PDF (created if absent)')
args = parser.parse_args()

event_num, file_num = args.event_num, args.file_num

import os
os.makedirs(args.output_dir, exist_ok=True)

IDEALEFF = False
filename = f'{args.output_dir}/event_display_{file_num}_{event_num:04d}.pdf'
# filename = f'event_displays/failing_5/event_display_{file_num}_{event_num:04d}.pdf'

def generate_cluster_colors(n):
    """Generate n perceptually distinct colors using golden-ratio HSV spacing."""
    import colorsys
    golden_ratio = 0.618033988749895
    hue = 0.0
    colors = []
    for _ in range(n):
        rgb = colorsys.hsv_to_rgb(hue, 0.75, 0.92)
        colors.append('#{:02x}{:02x}{:02x}'.format(
            int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255)
        ))
        hue = (hue + golden_ratio) % 1.0
    return colors

ANA_FILE = f'../ntuple-hgtd/user.mcardiff.45809429.Output._{file_num}.SuperNtuple.root'
tree = uproot.open(ANA_FILE)["ntuple"]

branch = tree.arrays([
    # vertex properties
    'RecoVtx_z'           ,
    'RecoVtx_time'        ,
    'RecoVtx_isHS'        ,
    'TruthVtx_z'          ,
    'TruthVtx_time'       ,
    'TruthVtx_isHS'       ,
    'RecoVtx_hasValidTime',
    # track propertries
    'Track_z0'          ,
    'Track_var_z0'      ,
    'Track_pt'          , 
    'Track_eta'         ,
    'Track_theta'       , 
    'Track_phi'         ,
    'Track_time'        , 
    'Track_timeRes'     ,
    'Track_qOverP'      , 
    'Track_quality'     ,
    'Track_nHGTDHits'   ,
    'Track_hasValidTime',
    'Track_truthVtx_idx',
    # reco jet properties
    'AntiKt4EMTopoJets_pt' ,
    'AntiKt4EMTopoJets_eta',
    'AntiKt4EMTopoJets_phi',
    'AntiKt4EMTopoJets_truthHSJet_idx',
])

reco_hs_z = branch.RecoVtx_z[event_num][0]
reco_hs_t = branch.RecoVtx_time[event_num][0]

truth_hs_z = branch.TruthVtx_z[event_num][0]
truth_hs_t = branch.TruthVtx_time[event_num][0]

# --- Track and Jet Processing ---
connected_tracks = []
for (idx, eta) in enumerate(branch.Track_eta[event_num]):
    in_hgtd = abs(eta) > 2.38 and abs(eta) < 4.0
    time_valid = True if IDEALEFF else branch.Track_hasValidTime[event_num][idx] == 1
    track_quality = branch.Track_quality[event_num][idx] == 1
    dz = branch.Track_z0[event_num][idx] - branch.RecoVtx_z[event_num][0]
    nsigma = abs(dz / np.sqrt(branch.Track_var_z0[event_num][idx]))
    if nsigma < 3.0 and branch.Track_pt[event_num][idx] > 1.0 and in_hgtd and time_valid:
        connected_tracks.append(idx)

track_info = []
for idx in connected_tracks:
    p         = abs(1 / branch.Track_qOverP[event_num][idx])
    track_eta = branch.Track_eta[event_num][idx]
    track_phi = branch.Track_phi[event_num][idx]
    track_pT  = branch.Track_pt[event_num][idx]
    track_z0  = branch.Track_z0[event_num][idx]
    track_vtx = branch.Track_truthVtx_idx[event_num][idx]
    status    = branch.TruthVtx_isHS[event_num][track_vtx]

    pz = track_pT * np.sinh(track_eta)
    signX = np.sign(track_eta) if track_eta != 0 else 1
    signY = np.sign(np.sin(track_phi)) if np.sin(track_phi) != 0 else 1
    theta = np.arctan(track_pT / abs(pz))
    x = (track_pT / 2) * np.cos(theta) * signX
    y = (track_pT / 2) * np.sin(theta) * signY
    track_info.append({'z0':track_z0, 'x':x, 'y':y, 'stat':status, 'idx':idx})

jet_info = []
for idx, jet_pt in enumerate(branch.AntiKt4EMTopoJets_pt[event_num]):
    if jet_pt < 30.0:
        continue
    jet_eta = branch.AntiKt4EMTopoJets_eta[event_num][idx]
    jet_phi = branch.AntiKt4EMTopoJets_phi[event_num][idx]
    isHS = len(branch.AntiKt4EMTopoJets_truthHSJet_idx[event_num][idx])

    pz = jet_pt * np.sinh(jet_eta)
    signX = np.sign(jet_eta) if jet_eta != 0 else 1
    signY = np.sign(np.sin(jet_phi)) if np.sin(jet_phi) != 0 else 1
    theta = np.arctan(jet_pt / abs(pz))
    x_jet = (jet_pt / 40) * np.cos(theta) * signX
    y_jet = (jet_pt / 40) * np.sin(theta) * signY
    jet_info.append({'pt':jet_pt, 'eta':jet_eta, 'phi':jet_phi, 'isHS':isHS,
                     'x':x_jet, 'y':y_jet, 'idx':idx})


# --- ROOT Macro Execution and Clustering ---
try:
    MACRO_CALL = f'runHGTD_Clustering.cxx("{file_num}",{event_num})'
    print(MACRO_CALL)
    result = subprocess.run(['root', '-l', '-q', '-b', MACRO_CALL],
                            check=True, capture_output=True, text=True)

    print(result.stdout)
    track_clusters = []
    cluster_times, cluster_zs = [], []
    cluster_trkptz_scores, cluster_waves_scores = [], []
    track_times = []
    current_block_idx = []
    current_block_times = []
    current_trkptz = None
    current_waves  = None

    for line in result.stdout.splitlines():
        line = line.strip()
        if line == "---------":
            if current_block_idx:
                track_clusters.append(current_block_idx)
                track_times.append(current_block_times)
                cluster_trkptz_scores.append(current_trkptz)
                cluster_waves_scores.append(current_waves)
                current_block_idx = []
                current_block_times = []
                current_trkptz = None
                current_waves  = None
            continue
        if line.startswith("t:"):
            cl_time = float(line[2:])
            if np.abs(cl_time - truth_hs_t) < 1000:
                cluster_times.append(cl_time)
            else:
                continue
        elif line.startswith("score_trkptz:"):
            current_trkptz = float(line.split(":", 1)[1])
        elif line.startswith("score_waves:"):
            current_waves = float(line.split(":", 1)[1])
        else:
            try:
                tup = line.split(",")
                trk_idx = int(tup[0])
                trk_time = float(tup[1])
                if np.abs(trk_time - truth_hs_t) < 1000:
                    current_block_idx.append(trk_idx)
                    current_block_times.append(trk_time)
            except (ValueError, IndexError):
                continue
    if current_block_idx:
        track_clusters.append(current_block_idx)
        track_times.append(current_block_times)
        cluster_trkptz_scores.append(current_trkptz)
        cluster_waves_scores.append(current_waves)
except subprocess.CalledProcessError as e:
    print(f"Error executing root script: {e.stderr}")
    track_clusters = []
    cluster_times = []

# --- Data for Histograms and Plotting ---
hist_times, time_errors = [], []
hs_times  , hs_zs       = [], []
hist_zs   , z_errors    = [], []
pt_wghts = []

for i_trk, cluster in enumerate(track_clusters):
    # Use list comprehensions to extract relevant data for the cluster
    this_zcluster = [branch.Track_z0[event_num][idx] for idx in cluster]
    this_var_z0 = [branch.Track_var_z0[event_num][idx] for idx in cluster]
    zbar_num = np.sum(np.array(this_zcluster) / np.array(this_var_z0))
    zbar_den = np.sum(1 / np.array(this_var_z0))

    hs_times.extend(track_times[i_trk][j_trk]
                    for (j_trk,idx) in enumerate(cluster)
                    if branch.Track_truthVtx_idx[event_num][idx] == 0)
    hs_zs.extend(branch.Track_z0[event_num][idx]
                 for idx in cluster
                 if branch.Track_truthVtx_idx[event_num][idx] == 0)

    hist_times.append(track_times[i_trk])
    hist_zs.append(this_zcluster)

    if IDEALEFF:
        time_errors.append([branch.Track_timeRes[event_num][idx]/np.sqrt(branch.Track_nHGTDHits[event_num][idx])
                            if branch.Track_hasValidTime[event_num][idx]==1 and branch.Track_nHGTDHits[event_num][idx] > 0 else 30
                            for idx in cluster])
    else:
        time_errors.append([branch.Track_timeRes[event_num][idx] for idx in cluster])
    z_errors.append(np.sqrt(this_var_z0))
    cluster_zs.append(zbar_num / zbar_den)
    pt_wghts.append([branch.Track_pt[event_num][idx] for idx in cluster])

print(hs_times)
print(hs_zs)

### NEW STUFF
pu_removed_tracks = []
for idx_trk in connected_tracks:
    z0_trk = branch.Track_z0[event_num][idx_trk]
    var_z0_trk = branch.Track_var_z0[event_num][idx_trk]

    closest_nsigma = np.inf
    closest_reco_vtx = -1
    for (idx_vtx, z_vtx) in enumerate(branch.RecoVtx_z[event_num]):
        nsigma = np.abs(z_vtx - z0_trk)/np.sqrt(var_z0_trk)
        if nsigma < closest_nsigma:
            closest_nsigma = nsigma
            closest_reco_vtx = idx_vtx
            
    if closest_reco_vtx != 0:
        pu_removed_tracks.append(idx_trk)

print(f'REMOVED TRACKS ', pu_removed_tracks)
### END NEW STUFF

# --- Plotting Functions and Generation ---
def draw_eta_reference_lines(ax, z_pos=reco_hs_z, eta_ref=2.4, line_length=50):
    """Draw Lines corresponding to |eta| = eta_ref wrt provided z position"""
    theta_ref = 2 * np.arctan(np.exp(-abs(eta_ref)))
    x_ref = line_length * np.cos(theta_ref)
    y_ref = line_length * np.sin(theta_ref)
    props = {'linestyle':'dotted', 'color':'lightgrey'}
    for dx in [x_ref, -x_ref]:
        for dy in [y_ref, -y_ref]:
            ax.plot([z_pos, z_pos + dx], [0, dy], **props)

def plot_rz_display(ax, track_info_list, jet_info_list):
    """Plot Event Display in R-Z Plane"""
    # Plot event display (z vs R)
    for track in track_info_list:
        ax.plot([track['z0'], track['z0'] + track['x']],[0, track['y']],
                color='blue' if track['stat'] == 1 else 'red',
                linestyle='solid')

    for (jet_i, jet_tup) in enumerate(jet_info_list):
        highlighted = (args.jet_idx is not None and jet_tup.get('idx') == args.jet_idx)
        jet_color = 'orange' if highlighted else ('green' if jet_tup['isHS'] >= 1 else 'grey')
        x_off1, y_off1 = jet_tup['x']-0.15*jet_tup['y'], jet_tup['y']+0.15*jet_tup['x']
        x_off2, y_off2 = jet_tup['x']+0.15*jet_tup['y'], jet_tup['y']-0.15*jet_tup['x']
        ax.fill([reco_hs_z, reco_hs_z + x_off1, reco_hs_z + x_off2], [0, y_off1, y_off2],
                color=jet_color, alpha=0.7 if highlighted else 0.5)

        txt_color = 'orange' if highlighted else ('green' if jet_tup['isHS'] >= 1 else 'black')
        label = f"Jet {jet_i+1}: $p_T$={jet_tup['pt']:.0f} GeV, $\eta$={jet_tup['eta']:.1f}"
        if highlighted:
            label += "  ← target"
        ax.text(reco_hs_z - 6.8, 0.9 - (1.2 + jet_i*0.1),
                label, weight='bold', fontsize=12, color=txt_color)

    draw_eta_reference_lines(ax, reco_hs_z, 2.4)
    ax.set_ylim(-1.0, 1.0)
    ax.set_xlim(reco_hs_z - 7.0, reco_hs_z + 7.0)
    ax.set_title(f'Event# {event_num}: Reco Vertex# {0}')
    ax.set_xlabel('Z [mm]')
    ax.set_yticks([])

    # Add vertices and labels
    ax.scatter(branch.RecoVtx_z[event_num], [0]*len(branch.RecoVtx_z[event_num]),
               marker='o', s=100,
               color=['blue' if z == reco_hs_z else 'black' for z in branch.RecoVtx_z[event_num]])
    ax.axhline(y=0.0, color='black', linestyle='--')

    ax.scatter(branch.TruthVtx_z[event_num], [-0.75]*len(branch.TruthVtx_z[event_num]),
               marker='|', s=100,
               color=['blue' if z == truth_hs_z else 'black' for z in branch.TruthVtx_z[event_num]])
    ax.axhline(y=-0.75, color='black', linestyle='--')

    ax.text(reco_hs_z + 4.5, 0.05, 'Reco vertices', fontsize=12)
    ax.text(reco_hs_z + 4.5, 0.05-0.75, 'Truth vertices', fontsize=12)
    ax.text(reco_hs_z + 2.0, 0.85, 'ATLAS Simulation Internal',
            fontsize=16, weight='bold', style='italic')

    # Target jet annotation (bottom-left, only when --jet_idx is given)
    if args.jet_idx is not None:
        has_rpt = args.rpt_hgtd is not None and args.rpt_mine is not None
        lines = [f'Target jet (idx={args.jet_idx})']
        if args.jet_label is not None:
            lines.append(f'  Truth identity: {args.jet_label}')
        if has_rpt:
            lines.append(f'  HGTD $R_{{pT}}$ = {args.rpt_hgtd:.3f}')
            lines.append(f'  My algo $R_{{pT}}$ = {args.rpt_mine:.3f}')
        for i, line in enumerate(lines):
            ax.text(reco_hs_z + 2.0, 0.72 - i * 0.13, line,
                    fontsize=11, color='orange',
                    weight='bold' if i == 0 else 'normal')

    # Add RZ display labels and legend
    ax.text(reco_hs_z-6.8, 0.9,
            f"Reco HS (z,t) = ({reco_hs_z:.1f} mm, {reco_hs_t:.1f} ps)",
            weight='bold', fontsize=12)
    ax.text(reco_hs_z-6.8, 0.9 - 0.1,
            f"Truth HS (z,t) = ({truth_hs_z:.1f} mm, {truth_hs_t:.1f} ps)",
            weight='bold', fontsize=12)

    # Jet and track legend
    legend_handles = [
        mlines.Line2D([], [], color='blue'),
        mlines.Line2D([], [], color='red'),
        mpatches.Rectangle((0, 0), 1, 1, color='green', alpha=0.5),
        mpatches.Rectangle((0, 0), 1, 1, color='grey',  alpha=0.5),
    ]
    legend_labels = ['Hard Scatter', 'Pile-Up', 'HS jet', 'PU jet']
    if args.jet_idx is not None:
        legend_handles.append(mpatches.Rectangle((0, 0), 1, 1, color='orange', alpha=0.7))
        legend_labels.append('Target jet')
    ax.legend(legend_handles, legend_labels,
              loc='upper left', title='Track and Jet Types',
              bbox_to_anchor=(0.0, 0.9))

def add_annotation(ax, truth_text_y, reco_text_y):
    """Add Arrows denoting the various important times"""
    y_min, y_max = ax.get_ylim()
    y_step = np.abs(y_max-y_min)/25
    x_min, x_max = ax.get_xlim()
    x_step = np.abs(x_max-x_min)/50
    ax.annotate("Truth HS Time",
                xytext=(truth_hs_t-5*x_step, y_min+8*y_step),
                xy=(truth_hs_t, truth_text_y), xycoords='data',
                ha='left', va='top', color='blue',
                arrowprops={'arrowstyle':'->','color':'blue','lw':2})
    ax.axvline(x=truth_hs_t, ymin=y_min, ymax=y_max, color='blue', lw=2)

    if branch.RecoVtx_hasValidTime[event_num][0] == 1:
        ax.annotate("HGTD Time", xytext=(reco_hs_t-5*x_step, y_min-2*y_step),
                    xy=(reco_hs_t, reco_text_y), xycoords='data',
                    ha='left', va='top', color='black',
                    arrowprops={'arrowstyle':'->','color':'black','lw':2})
        ax.axvline(x=reco_hs_t, ymin=y_min, ymax=y_max, color='black', lw=2)

    if args.extra_time:
        ax.annotate("My Time", xytext=(args.extra_time-2.5*x_step, y_min-6*y_step),
                    xy=(args.extra_time, reco_text_y),
                    ha='left', va='top', color='black',
                    arrowprops={'arrowstyle':'->','color':'black','lw':2})

cluster_colors = generate_cluster_colors(len(track_clusters))

random.seed(42069)
calo_time = random.gauss(truth_hs_t,90)
calo_min, calo_max = calo_time-3*90, calo_time+3*90
calo_label = f'Calo Time Exclusion\n[{calo_min:.2f},{calo_max:.2f}]'

all_times = [reco_hs_t, truth_hs_t] + list(chain.from_iterable(hist_times))
extended_min_time = min(all_times) - 0.05 * (max(all_times) - min(all_times))
extended_max_time = max(all_times) + 0.05 * (max(all_times) - min(all_times))

# Rank clusters by TRKPTZ score (parsed from C++ output); fall back to
# Python-recomputed sum(pT)*exp(-1.5|Δz|) if the parsed score is missing.
def _trkptz_or_fallback(i):
    s = cluster_trkptz_scores[i] if i < len(cluster_trkptz_scores) else None
    if s is not None:
        return s
    return sum(pt_wghts[i]) * np.exp(-1.5 * np.abs(reco_hs_z - cluster_zs[i]))

# Top 4 by TRKPTZ score (descending)
ranked_indices = sorted(range(len(pt_wghts)), key=_trkptz_or_fallback, reverse=True)
legend_indices = ranked_indices[:4]

# Always include the cluster whose time is closest to the truth HS vertex time
truth_closest_idx = min(range(len(cluster_times)), key=lambda i: np.abs(cluster_times[i] - truth_hs_t))
if truth_closest_idx not in legend_indices:
    legend_indices.append(truth_closest_idx)

legend_index_set = set(legend_indices)

weight_legend_patches = []
for i_cluster, weights in enumerate(pt_wghts):
    if i_cluster not in legend_index_set:
        continue
    t_label = f't = {cluster_times[i_cluster]:.2f} ps'
    suffix  = '  [ΔT closest]' if i_cluster == truth_closest_idx and truth_closest_idx not in ranked_indices[:4] else ''
    trkptz_s = cluster_trkptz_scores[i_cluster] if i_cluster < len(cluster_trkptz_scores) else None
    waves_s  = cluster_waves_scores[i_cluster]  if i_cluster < len(cluster_waves_scores)  else None
    if trkptz_s is not None and waves_s is not None:
        score_label = f'TRKPTZ = {trkptz_s:.3g}  |  WAVeS = {waves_s:.3g}'
    elif trkptz_s is not None:
        score_label = f'TRKPTZ = {trkptz_s:.3g}'
    else:
        dz = np.abs(reco_hs_z - cluster_zs[i_cluster])
        score_label = f'score = {sum(weights) * np.exp(-1.5 * dz):.3f}'
    FULL_LABEL = f'{score_label}\n{t_label}{suffix}'
    patch = mpatches.Patch(facecolor=cluster_colors[i_cluster], label=FULL_LABEL)
    weight_legend_patches.append(patch)

with PdfPages(filename) as pdf:
    # --- Page 1: R-Z and T-Z Track Dist ---
    fig1, axs1 = plt.subplots(2, 1, figsize=(12, 12))
    event_display1, time_histo1 = axs1

    # ---- Upper Region: R-Z Display ----
    plot_rz_display(event_display1, track_info, jet_info)

    # ---- Lower Region: Z-T Display ----
    full_info = zip(hist_times, time_errors, hist_zs, z_errors, cluster_colors)
    for (times, errtime, zs, errz, color) in full_info:
        time_histo1.errorbar(
            times, zs, fmt='.',
            xerr=errtime, yerr=errz,
            lw=1.5, alpha=1.0, color=color,
            zorder=3)

    # Hard Scatter times
    time_histo1.scatter(hs_times, hs_zs, marker='o', edgecolors='black',
                        s=50, color='blue', label='Hard Scatter Tracks', zorder=15)

    # Add Reco Vertex Z
    time_histo1.plot(list(time_histo1.get_xlim()), [reco_hs_z]*2,
                     label='Reco Vertex Z', color='black', lw=2, linestyle='--')

    # Cluster Times
    time_histo1.scatter(cluster_times, cluster_zs,
                        marker="*", edgecolors='black',
                        s=500, color=cluster_colors,
                        label='Cluster Positions', zorder=2)

    # Time Annotations
    add_annotation(time_histo1, truth_hs_z, reco_hs_z)

    # Calo Time Exclusion Region
    time_histo1.axvspan(time_histo1.get_xlim()[0], calo_min,
                        color='grey', alpha=0.2, zorder=1,
                        label=calo_label)

    time_histo1.axvspan(time_histo1.get_xlim()[1], calo_max,
                        color='grey', alpha=0.2, zorder=1)

    # Set Limits and Labels
    time_histo1.set_xlim(extended_min_time, extended_max_time)
    time_histo1.set_xlabel('Time (ps)')
    time_histo1.set_ylabel('Z (mm)')
    time_histo1.set_title('Z-T Event Display')

    # Legends
    cluster_legend = time_histo1.legend(
        handles=weight_legend_patches,
        loc='upper right',
        title="Cluster Weights",
        fontsize='small',
        title_fontsize='small')

    time_histo1.add_artist(cluster_legend)
    time_histo1.legend(loc='upper left')

    pdf.savefig(fig1)
    plt.close(fig1)

    # --- Page 2: RZ and Time Histogram ---
    fig2, axs2 = plt.subplots(2, 1, figsize=(12, 12))
    event_display2, time_histo2 = axs2

    # ---- Upper Region: R-Z Display ----
    plot_rz_display(event_display2, track_info, jet_info)

    # ---- Lower Region: Z-T Display ----
    histvals, bin_edges, patches = time_histo2.hist(
        hist_times, bins=50, stacked=True,
        range=(extended_min_time, extended_max_time),
        weights=pt_wghts, label='Track Time',
        lw=1.5, alpha=1.0, color=cluster_colors)

    max_val = histvals.max()
    time_histo2.set_ylim(0, 1.1 * max_val)

    flat_times = list(chain.from_iterable(hist_times))
    flat_pt = list(chain.from_iterable(pt_wghts))

    # Hard Scatter times
    for t in hs_times:
        time_histo2.text(t, 0, '/', ha='center', va='top', fontsize=20, color='blue')

    # Hatch the bin for the highest pT track
    max_pt_idx = np.argmax(flat_pt)
    bin_idx = np.digitize(flat_times[max_pt_idx], bin_edges) - 1

    max_pt_rect = mpatches.Rectangle(
        (bin_edges[bin_idx], 0),
        bin_edges[bin_idx+1] - bin_edges[bin_idx],
        max(flat_pt), # Fixed height as a fraction of max bin height
        linewidth=1.5,
        facecolor='none',
        edgecolor='black',
        hatch='///',
        label='Highest pT Track',
        zorder=3)
    time_histo2.add_patch(max_pt_rect)

    # Cluster Times
    time_histo2.scatter(cluster_times, [0.1*max_val]*len(cluster_times),
                        marker="*", edgecolors='black',
                        s=500, color=cluster_colors,
                        label='Cluster Positions', zorder=4)

    # Time Annotations
    add_annotation(time_histo2, 0.1*max_val, 0.1*max_val)

    # Calo Time Exclusion Region
    time_histo2.axvspan(time_histo2.get_xlim()[0], calo_min,
                        color='grey', alpha=0.2, zorder=1,
                        label=calo_label)

    time_histo2.axvspan(time_histo2.get_xlim()[1], calo_max,
                        color='grey', alpha=0.2, zorder=1)

    # set labels and such
    time_histo2.set_xlim(extended_min_time, extended_max_time)
    time_histo2.set_xlabel('Time (ps)')
    time_histo2.set_ylabel('Track $p_T$')
    time_histo2.set_title('Time Histogram')

    # Legends
    cluster_legend = time_histo2.legend(
        handles=weight_legend_patches,
        loc='upper right',
        title="Cluster Weights",
        fontsize='small',
        title_fontsize='small')

    time_histo2.add_artist(cluster_legend)
    time_histo2.legend(loc='upper left')

    pdf.savefig(fig2)
    plt.close(fig2)

    # --- Page 3: η-φ cluster view (split: negative η left, positive η right) ---
    ETA_MIN, ETA_MAX = 2.38, 4.0
    fig3, (ax_neg, ax_pos) = plt.subplots(1, 2, figsize=(14, 8), sharey=True)
    fig3.subplots_adjust(bottom=0.22, wspace=0.05)

    def _eta_ax(eta):
        return ax_neg if eta < 0 else ax_pos

    # Per-cluster tracks in η-φ (signed eta, correct panel)
    for i_cl, cluster in enumerate(track_clusters):
        col = cluster_colors[i_cl]
        for idx in cluster:
            eta = float(branch.Track_eta[event_num][idx])
            phi = float(branch.Track_phi[event_num][idx])
            pt  = float(branch.Track_pt[event_num][idx])
            is_hs = branch.Track_truthVtx_idx[event_num][idx] == 0
            _eta_ax(eta).scatter(eta, phi,
                        s=max(20, pt * 8),
                        color=col,
                        edgecolors='blue' if is_hs else 'none',
                        linewidths=0.8 if is_hs else 0,
                        zorder=4 if is_hs else 3,
                        marker='o')

    # Cluster centroids: signed mean η and circular mean φ
    for i_cl, cluster in enumerate(track_clusters):
        col  = cluster_colors[i_cl]
        etas = np.array([float(branch.Track_eta[event_num][idx]) for idx in cluster])
        phis = np.array([float(branch.Track_phi[event_num][idx]) for idx in cluster])
        eta_c = np.mean(etas)
        phi_c = np.arctan2(np.mean(np.sin(phis)), np.mean(np.cos(phis)))
        _eta_ax(eta_c).scatter(eta_c, phi_c,
                    s=350, marker='*', color=col,
                    edgecolors='black', linewidths=0.7, zorder=6)

    # Jets as ΔR = 0.4 circles (placed in the correct signed-η panel)
    for jet in jet_info:
        jet_eta = jet['eta']
        if abs(jet_eta) + 0.4 < ETA_MIN or abs(jet_eta) - 0.4 > ETA_MAX:
            continue
        highlighted = (args.jet_idx is not None and jet.get('idx') == args.jet_idx)
        jet_color = 'orange' if highlighted else ('green' if jet['isHS'] else 'grey')
        lw = 3.0 if highlighted else 1.5
        ax_j = _eta_ax(jet_eta)
        circle = plt.Circle((jet_eta, jet['phi']), 0.4,
                             color=jet_color, fill=False,
                             linewidth=lw, linestyle='-',
                             alpha=0.9 if highlighted else 0.7, zorder=5)
        ax_j.add_patch(circle)
        label = f" {jet['pt']:.0f} GeV" + (" ← target" if highlighted else "")
        ax_j.text(jet_eta + 0.02 * np.sign(jet_eta), jet['phi'],
                  label, fontsize=8, color=jet_color, va='center', zorder=6)

    # Axis limits and labels
    ax_neg.set_xlim(-ETA_MAX, -ETA_MIN)
    ax_pos.set_xlim(ETA_MIN, ETA_MAX)
    for ax in (ax_neg, ax_pos):
        ax.set_ylim(-np.pi, np.pi)
        ax.set_yticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
        ax.set_yticklabels(['-π', '-π/2', '0', 'π/2', 'π'])
        ax.set_xlabel('η', fontsize=13)
        ax.grid(True, alpha=0.3)
    ax_neg.set_ylabel('φ  (rad)', fontsize=13)
    ax_neg.set_title('Negative η  (backward)', fontsize=12)
    ax_pos.set_title('Positive η  (forward)', fontsize=12)
    fig3.suptitle(f'Event# {event_num}: η-φ Cluster View  (HGTD acceptance)', fontsize=13, y=1.01)
    ax_neg.text(0.02, 0.98, 'ATLAS Simulation Internal',
                transform=ax_neg.transAxes, fontsize=11,
                weight='bold', style='italic', va='top')

    # Legend below axes
    leg_handles = [
        mlines.Line2D([], [], marker='o', color='none',
                      markerfacecolor='grey', markeredgecolor='blue',
                      markersize=8, label='Hard Scatter track (blue edge)'),
        mlines.Line2D([], [], marker='o', color='none',
                      markerfacecolor='grey', markeredgecolor='none',
                      markersize=8, label='Pile-up track'),
        mlines.Line2D([], [], marker='*', color='none',
                      markerfacecolor='grey', markeredgecolor='black',
                      markersize=12, label='Cluster centroid (mean η, φ)'),
        mpatches.Patch(facecolor='green',  alpha=0.6, label='HS jet (ΔR=0.4)'),
        mpatches.Patch(facecolor='grey',   alpha=0.5, label='PU jet (ΔR=0.4)'),
    ]
    if args.jet_idx is not None:
        leg_handles.append(mpatches.Patch(facecolor='none', edgecolor='orange',
                                          linewidth=2.5, label='Target jet (ΔR=0.4)'))
    fig3.legend(handles=leg_handles,
                loc='upper center', bbox_to_anchor=(0.5, -0.02),
                ncol=5, fontsize=9, framealpha=0.9)

    pdf.savefig(fig3, bbox_inches='tight')
    plt.close(fig3)
