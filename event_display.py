"""
EVENT DISPLAY CODE
CREATED BY: wasikul.islam@cern.ch
MODIFIED FOR USE BY: mcardiff@brandeis.edu
"""

import argparse
import random
import subprocess
from operator import __sub__
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
args = parser.parse_args()

event_num, file_num = args.event_num, args.file_num

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
    time_valid = branch.Track_hasValidTime[event_num][idx] == 1
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
    theta = np.atan(track_pT / abs(pz))
    x = (track_pT / 2) * np.cos(theta) * signX
    y = (track_pT / 2) * np.sin(theta) * signY
    track_info.append({'z0':track_z0, 'x':x, 'y':y, 'stat':status})

jet_info = []
for jet_pt in branch.AntiKt4EMTopoJets_pt[event_num]:
    idx = branch.AntiKt4EMTopoJets_pt[event_num].tolist().index(jet_pt)
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
    jet_info.append({'pt':jet_pt, 'eta':jet_eta, 'phi':jet_phi, 'isHS':isHS, 'x':x_jet, 'y':y_jet})

# --- ROOT Macro Execution and Clustering ---
try:
    MACRO_CALL = f'runHGTD_Clustering.cxx("{file_num}",{event_num})'
    result = subprocess.run(['root', '-l', '-q', '-b', MACRO_CALL],
                            check=True, capture_output=True, text=True)
    track_clusters = []
    cluster_times, cluster_zs = [], []
    track_times = []
    current_block_idx = []
    current_block_times = []

    for line in result.stdout.splitlines():
        line = line.strip()
        if line == "---------":
            if current_block_idx:
                track_clusters.append(current_block_idx)
                track_times.append(current_block_times)
                current_block_idx = []
                current_block_times = []
            continue
        if "t:" in line:
            cluster_times.append(float(line[2:]))
        try:
            tup = line.split(",")
            trk_idx = int(tup[0])
            trk_time = float(tup[1])
            current_block_idx.append(trk_idx)
            current_block_times.append(trk_time)
        except (ValueError, IndexError):
            continue
    if current_block_idx:
        track_clusters.append(current_block_idx)
        track_times.append(current_block_times)
except subprocess.CalledProcessError as e:
    print(f"Error executing root script: {e.stderr}")
    track_clusters = []
    cluster_times = []

# --- Data for Histograms and Plotting ---
hist_times, time_errors = [], []
hist_zs   , z_errors    = [], []
hs_times  , hs_zs       = [], []
pt_wghts = []

for i_trk, cluster in enumerate(track_clusters):
    # Use list comprehensions to extract relevant data for the cluster
    this_zcluster = [branch.Track_z0[event_num][idx] for idx in cluster]
    this_var_z0 = [branch.Track_var_z0[event_num][idx] for idx in cluster]
    zbar_num = np.sum(np.array(this_zcluster) / np.array(this_var_z0))
    zbar_den = np.sum(1 / np.array(this_var_z0))

    hist_times.append(track_times[i_trk])
    hist_zs.append(this_zcluster)
    time_errors.append([branch.Track_timeRes[event_num][idx] for idx in cluster])
    z_errors.append(np.sqrt(this_var_z0))
    cluster_zs.append(zbar_num / zbar_den)
    pt_wghts.append([branch.Track_pt[event_num][idx] for idx in cluster])

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
                color='blue' if track['stat'] == 1 else 'red')

    for (jet_i, jet_tup) in enumerate(jet_info_list):
        jet_color = 'green' if jet_tup['isHS'] >= 1 else 'grey'
        x_off1, y_off1 = jet_tup['x']-0.15*jet_tup['y'], jet_tup['y']+0.15*jet_tup['x']
        x_off2, y_off2 = jet_tup['x']+0.15*jet_tup['y'], jet_tup['y']-0.15*jet_tup['x']
        ax.fill([reco_hs_z, reco_hs_z + x_off1, reco_hs_z + x_off2], [0, y_off1, y_off2],
                color=jet_color, alpha=0.5)

        txt_color = 'green' if jet_tup['isHS'] >= 1 else 'black'
        ax.text(reco_hs_z - 6.8, 0.9 - (1.2 + jet_i*0.1),
                f"Jet {jet_i+1}: $p_T$={jet_tup['pt']:.0f} GeV, $\eta$={jet_tup['eta'] :.1f}",
                weight='bold', fontsize=12, color=txt_color)

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

    # Add RZ display labels and legend
    ax.text(reco_hs_z-6.8, 0.9,
            f"Reco HS (z,t) = ({reco_hs_z:.1f} mm, {reco_hs_t:.1f} ps)",
            weight='bold', fontsize=12)
    ax.text(reco_hs_z-6.8, 0.9 - 0.1,
            f"Truth HS (z,t) = ({truth_hs_z:.1f} mm, {truth_hs_t:.1f} ps)",
            weight='bold', fontsize=12)

    # Jet and track legend
    ax.legend([mlines.Line2D([], [], color='blue'),
               mlines.Line2D([], [], color='red'),
               mpatches.Rectangle((0, 0), 1, 1, color='green', alpha=0.5),
               mpatches.Rectangle((0, 0), 1, 1, color='grey', alpha=0.5)],
              ['Hard Scatter', 'Pile-Up', 'HS jet', 'PU jet'],
              loc='upper left', title='Track and Jet Types',
              bbox_to_anchor=(0.0, 0.9))

def add_annotation(ax, truth_text_y, reco_text_y):
    """Add Arrows denoting the various important times"""
    y_min, y_max = ax.get_ylim()
    y_step = np.abs(y_max-y_min)/25
    x_min, x_max = ax.get_xlim()
    x_step = np.abs(x_max-x_min)/50
    ax.annotate("Truth HS Time",
                xytext=(truth_hs_t-2.5*x_step, y_min-2*y_step),
                xy=(truth_hs_t, truth_text_y), xycoords='data',
                ha='left', va='top', color='blue',
                arrowprops={'arrowstyle':'->','color':'blue','lw':2})

    if branch.RecoVtx_hasValidTime[event_num][0] == 1:
        ax.annotate("Reco HS Time", xytext=(reco_hs_t-2.5*x_step, y_min-4*y_step),
                    xy=(reco_hs_t, reco_text_y), xycoords='data',
                    ha='left', va='top', color='black',
                    arrowprops={'arrowstyle':'->','color':'black','lw':2})

    if args.extra_time:
        ax.annotate("My Time", xytext=(args.extra_time-2.5*x_step, y_min-6*y_step),
                    xy=(args.extra_time, reco_text_y),
                    ha='left', va='top', color='black',
                    arrowprops={'arrowstyle':'->','color':'black','lw':2})

filename = f'./figs/nohs/event_display_{file_num}_{event_num:04d}.pdf'
colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628", "#f781bf", "#999999"]
cluster_colors = colors[:len(track_clusters)]

random.seed(42069)
calo_time = random.gauss(truth_hs_t,90)
calo_min, calo_max = calo_time-3*90, calo_time+3*90
calo_label = f'Calo Time Exclusion\n[{calo_min:.2f},{calo_max:.2f}]'

all_times = [reco_hs_t, truth_hs_t] + list(chain.from_iterable(hist_times))
extended_min_time = min(all_times) - 0.05 * (max(all_times) - min(all_times))
extended_max_time = max(all_times) + 0.05 * (max(all_times) - min(all_times))

weight_legend_patches = []
for i_cluster, weights in enumerate(pt_wghts):
    total_weight = sum(weights)
    pt_label = f'Σpt = {total_weight:.2f}'
    dz_label = f'dz={np.abs(reco_hs_z-cluster_zs[i_cluster]):.2f}'
    t_label  = f't={cluster_times[i_cluster]:.2f}'
    FULL_LABEL = f'{pt_label}\n{dz_label}\n{t_label}'
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
                       s=50, color='blue', label='Hard Scatter Tracks', zorder=10)

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
        time_histo2.text(t, 0, '/', ha='center', va='top',
                         fontsize=20, color='blue')

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
