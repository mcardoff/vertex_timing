####################################################################
######## Event display code for vertices ###########################
######## Contact for comments to wasikul.islam@cern.ch #############
####################################################################

import uproot
import matplotlib.pyplot as plt
import numpy as np
import math
from math import log, tan, cosh, sqrt, pi
from matplotlib.lines import Line2D  # Explicit import of Line2D
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.patheffects as path_effects
from matplotlib.backends.backend_pdf import PdfPages
import random
import subprocess
from itertools import chain

#############################################################
import argparse
parser = argparse.ArgumentParser(description='Process event and vertex data from ROOT file.')
parser.add_argument('--event_num', type=int, required=True, help='Event number to process')
parser.add_argument('--vtxID', type=int, required=True, help='Vertex ID to process', default=0)
parser.add_argument('--file_num', type=str, required=True, help='Specific ID of Filen to process')
parser.add_argument('--extra_time', type=float, required=False, help='Optional other selected time')
args = parser.parse_args()

event_num = args.event_num
vtxID = args.vtxID
file_num = args.file_num
#############################################################

root_file = uproot.open(f'../ntuple-hgtd/user.mcardiff.45809429.Output._{file_num}.SuperNtuple.root')
tree = root_file["ntuple"]

my_branches = tree.arrays([
    'TruthVtx_z',
    'TruthVtx_time',
    'TruthVtx_isHS',
    'RecoVtx_isHS',
    'RecoVtx_z',
    'RecoVtx_time',
    'RecoVtx_hasValidTime',
    'Track_truthVtx_idx',
    'RecoVtx_sumPt2',
    # 'RecoVtx_isPU',
    'RecoVtx_track_idx',
    'AntiKt4EMTopoJets_ghostTrack_idx',
    'Track_z0',
    'Track_pt',
    'Track_eta',
    'Track_time',
    'Track_timeRes',
    'Track_hasValidTime',
    'Track_qOverP',
    'Track_theta',
    'Track_phi',
    'Track_var_z0',
    'Track_charge',
    'AntiKt4EMTopoJets_pt',
    'AntiKt4EMTopoJets_eta',
    'AntiKt4EMTopoJets_phi',
    'AntiKt4EMTopoJets_truthHSJet_idx',
    'TruthHSJet_pt',
    'TruthHSJet_eta',
    'TruthHSJet_phi',
    "Track_quality",
    'Track_nearestVtx_idx',
    'Track_nearestVtx_sig',
    'Track_truthPart_idx',
    'TruthPart_prodVtx_time',
    'TruthPart_pt',
    'TruthPart_eta',
    'TruthPart_phi',
    'TruthPart_charge',
])

# my_branches = tree.arrays(['TruthVtx_z', 'RecoVtx_sumPt2', 'RecoVtx_z'])

vtx_z = my_branches.RecoVtx_z[event_num][vtxID]
reco_hs_z = my_branches.RecoVtx_z[event_num][0]

truth_z = my_branches.TruthVtx_z[event_num][vtxID]
truth_hs_z = my_branches.TruthVtx_z[event_num][0]

vtx_t = my_branches.RecoVtx_time[event_num][vtxID]
reco_hs_t = my_branches.RecoVtx_time[event_num][0]

truth_t = my_branches.TruthVtx_time[event_num][vtxID]
truth_hs_t = my_branches.TruthVtx_time[event_num][0]

sumpt = my_branches.RecoVtx_sumPt2[event_num][vtxID]

selected_HS_vtx_id = None
my_track_z0=[]


# Iterate through the vertex z-values and sumpt values
max_sumpt = -1
for idx, z in enumerate(my_branches.RecoVtx_z[event_num]):
    # print(f"Reco z[{idx}]: {z:.2f}")
    if abs(z - truth_z) <= 5:
        sumpt_value = my_branches.RecoVtx_sumPt2[event_num][idx]
        #print("reco vertex within 5 mm of truth z : ", z, "vtx id:", idx, "sumpt :", sumpt_value, "|(z - truth_z)| :", abs(z - truth_z))
        if  sumpt_value > max_sumpt:
            selected_HS_vtx_id = idx
            max_sumpt = my_branches.RecoVtx_sumPt2[event_num][idx]
            
print(f"python event_display_VBF_R25.py --file_num {file_num} --event_num {event_num} --vtxID {selected_HS_vtx_id}")


# Print the selected vertex and its sumpt value
if selected_HS_vtx_id is not None:
    print("Evt# ", event_num, f" Selected Vertex ID: {selected_HS_vtx_id}")
    #print(f"Sumpt of Selected Vertex: {max_sumpt}")
else:
    print("No vertices found within the specified range.")

if (vtxID == selected_HS_vtx_id):
    print("This is HS vertex")
else :
    print("This is PU vertex")

# Get the vertex z-coordinate
print("isHS =",my_branches.RecoVtx_isHS[event_num][vtxID],
      "RecoVtx Z = ", my_branches.RecoVtx_z[event_num][vtxID]
      # , " isPU =",my_branches.RecoVtx_isPU[event_num][vtxID]
      )
print("isTrueHS =",my_branches.TruthVtx_isHS[event_num][vtxID],
      "TruthVtx Z = ", my_branches.TruthVtx_z[event_num][vtxID]
      # , " isPU =",my_branches.RecoVtx_isPU[event_num][vtxID]
      )


closest_vertices = []
reco_vertices = []
for idx, z in enumerate(my_branches.RecoVtx_z[event_num]):
    if abs(z - vtx_z) <= 5:
        reco_vertices.append((z, idx))
        diff = abs(z - vtx_z)
        closest_vertices.append((idx, z, diff))

closest_vertices.sort(key=lambda x: x[2])
# Display the closest vertices
for idx, z, diff in closest_vertices:
    #print(f"Vertex ID: {idx}, Z: {z}, Difference from vtx_z: {diff}")
    if (idx==selected_HS_vtx_id):
        print("closest HS vertex# :", idx)
    
truth_vertices = []
for idx, z in enumerate(my_branches.TruthVtx_z[event_num]):
    # print(f"Truth z[{idx}]: {z:.2f}")
    if abs(z - vtx_z) <= 5:
        truth_vertices.append((z, idx))

# min_t_diff = 1e6
# min_idx = 0
# matched_time = 0
# for idx, t in enumerate(my_branches.TruthVtx_time[event_num]):
#     if (abs(t-reco_hs_t)) < min_t_diff:
#         min_t_diff = abs(t-reco_hs_t)
#         min_idx = idx
#         matched_time = t
# matched_z = my_branches.TruthVtx_z[event_num][min_idx]
# print(f"Min diff is with truthVtx[{min_idx}]\nReco HS (z,t):({reco_hs_z:.2f},{reco_hs_t:.2f})\nTruth (z,t):({matched_z:.2f},{matched_time:.2f})")
        
#####################################################

# Get all tracks within 3 sigma of selected vertex
# connected_tracks = my_branches.RecoVtx_track_idx[event_num][vtxID]
connected_tracks = []
for idx in range(len(my_branches.Track_z0[event_num])):
    inHGTD = abs(my_branches.Track_eta[event_num][idx]) > 2.38 and abs(my_branches.Track_eta[event_num][idx]) < 4.0
    hasTime = my_branches.Track_hasValidTime[event_num][idx] == 1
    dz = my_branches.Track_z0[event_num][idx] - my_branches.RecoVtx_z[event_num][vtxID]
    nsigma = abs(dz / sqrt(my_branches.Track_var_z0[event_num][idx]))
    isHS = my_branches.Track_truthVtx_idx[event_num][idx] == 0
    if(abs(nsigma) < 3.0 and my_branches.Track_pt[event_num][idx] > 1.0 and inHGTD and hasTime):
        connected_tracks.append(idx)

track_info = []
jet_info = []
num_jets = 0
new_sumpt=0

num_HS_tracks=0

for idx in connected_tracks:
    track_z0 = my_branches.Track_z0[event_num][idx]
    p = abs(1 / (my_branches.Track_qOverP[event_num][idx]))
    track_eta = -np.log(math.tan((my_branches.Track_theta[event_num][idx]) / 2))
    track_pT = (p / (np.cosh(track_eta))) / 1000
    track_phi = my_branches.Track_phi[event_num][idx]
    z0 = track_z0 - vtx_z
        
    my_track_z0.append(track_z0)
    
    # Loop over jets
    for j in range(len(my_branches.AntiKt4EMTopoJets_ghostTrack_idx[event_num])):
        if my_branches.AntiKt4EMTopoJets_pt[event_num][j] < 30.0:
            continue
        # njets_b1=njets_b1+1    
        trackPT0 = 0
        trackPT = 0
        jet_pt_05_rpt = 0
        jet_pt_01_rpt = 0

        for q in range(len(my_branches.AntiKt4EMTopoJets_ghostTrack_idx[event_num][j])):
            # Track selection cuts
            jdx = my_branches.AntiKt4EMTopoJets_ghostTrack_idx[event_num][j][q]
            p2 = abs(1 / my_branches.Track_qOverP[event_num][jdx])
            eta2 = -log(tan(my_branches.Track_theta[event_num][jdx] / 2))
            track_pT2 = p2 / cosh(eta2)
            pt2 = track_pT2 / 1000
    
            delz = my_branches.Track_z0[event_num][jdx] - my_branches.RecoVtx_z[event_num][vtxID]
            signi_cut = delz / sqrt(my_branches.Track_var_z0[event_num][jdx])
    
            if abs(signi_cut) > 3.0:
                # print("CONTINUING FROM INSIDE LOOP")
                continue
            trackPT += pt2
    
        Rpt = trackPT / my_branches.AntiKt4EMTopoJets_pt[event_num][j]
            
        jet_pt = my_branches.AntiKt4EMTopoJets_pt[event_num][j]
        jet_eta = my_branches.AntiKt4EMTopoJets_eta[event_num][j]
        jet_phi = my_branches.AntiKt4EMTopoJets_phi[event_num][j]
        jet_isHS = my_branches.AntiKt4EMTopoJets_truthHSJet_idx[event_num][j]
        jet_isHS_size = len(jet_isHS)
        jet_pz = jet_pt * np.sinh(jet_eta)

        jet_signX = np.sign(jet_eta) if jet_eta != 0 else 1  # Avoid division by zero
        jet_signY = np.sign(np.sin(jet_phi)) if np.sin(jet_phi) != 0 else 1

        jet_theta = np.arctan(jet_pt / abs(jet_pz))
        jet_x = (jet_pt / 40) * np.cos(jet_theta) * jet_signX
        jet_y = (jet_pt / 40) * np.sin(jet_theta) * jet_signY
        #print("Jet : ", " pt :", jet_pt, "eta :", jet_eta, "isHS :", jet_isHS_size)
        
        jet_tuple = (jet_pt, jet_eta, jet_phi, jet_isHS_size, jet_x, jet_y, Rpt)

        if jet_tuple not in jet_info:
            jet_info.append(jet_tuple)
        
        num_jets += 1  # Increment counter
    
    pz = track_pT * math.sinh(track_eta)
    signX = track_eta / abs(track_eta)
    signY = math.sin(track_phi) / abs(math.sin(track_phi))
    theta = math.atan(track_pT / abs(pz))
    x = (track_pT / 2) * math.cos(theta) * signX
    y = (track_pT / 2) * math.sin(theta) * signY

    Track_truthVtx_id = my_branches.Track_truthVtx_idx[event_num][idx]

    status = my_branches.TruthVtx_isHS[event_num][Track_truthVtx_id]
    dz_track = my_branches.Track_z0[event_num][idx] - vtx_z
    nsigma_track = dz_track / sqrt(my_branches.Track_var_z0[event_num][idx])

    isTruthHS = my_branches.TruthVtx_isHS[event_num][Track_truthVtx_id]==1
    hasTime = my_branches.Track_hasValidTime[event_num][idx] == 1

    # if (isTruthHS and (not hasTime)):
    #     status = 2
    # if ((not isTruthHS) and (not hasTime)):
    #     status = 3

    # if ((not hasTime) and isTruthHS):
    #     status = 4

    # if ((not hasTime) and (not isTruthHS)):
    #     status = 5
    

    truthpart_idx = my_branches.Track_truthPart_idx[event_num][idx]
    truthpart_pt = my_branches.TruthPart_pt[event_num][truthpart_idx]
    truthpart_eta = my_branches.TruthPart_eta[event_num][truthpart_idx]
    
    # goodMatch = (abs(track_pT - truthpart_pt)/truthpart_pt < 0.1) and (abs(track_eta - truthpart_eta) < 0.01) and (my_branches.Track_charge[event_num][idx] == my_branches.TruthPart_charge[event_num][truthpart_idx]);
    # if (not goodMatch):
    #     print("FAKE?")
    #     status = 4
    near = my_branches.Track_nearestVtx_idx[event_num][idx]
    nsigma_near = np.abs(my_branches.Track_nearestVtx_sig[event_num][idx])
    nearcut = not (near != vtxID and nsigma_near < 3.0)

    # if ((not nearcut) and isTruthHS):
    #     status = 2
    # elif ((not nearcut) and (not isTruthHS)):
    #     status = 3
    
    # if (nsigma_track < 3.0 and isTruthHS):
    if (np.abs(nsigma_track) < 3.0):
        track_info.append([vtx_z, z0, x, y, status])
        num_HS_tracks=num_HS_tracks+1

print("number of HS tracks : ", num_HS_tracks)

# print("max z0: ", max(track_info, key=lambda t: t[0]+t[1])[1])
# print("min z0: ", min(track_info, key=lambda t: t[0]+t[1])[1])

random.seed(42069)
guess_time = random.gauss(truth_hs_t,90)
# ONLY USE FOR SPECIFIC EVENT
# track_clusters = [[1384, 894, 994, 1236, 1511, 1645, 1717, 1723, 1739 ],
                  # [1012, 817, 1470, 1686, 1722],
                  # [1482, 1251]]
# cluster_times = [-29.2673, 137.551, 300.575]
track_clusters = []
track_times = []
cluster_times = []
use_smeared_times = False # True
try:
    macro_call = f'runHGTD_Clustering.cxx("{file_num}",{event_num},{"true" if use_smeared_times else "false"},true)'
    result = subprocess.run(['root', '-l', '-q', '-b', macro_call],
                            check=True,
                            capture_output=True,
                            text=True)
    current_block_idx = []
    current_block_times = []
    print('STDOUT FROM MACRO: ')
    print(result.stdout)

    for line in result.stdout.splitlines():
        line = line.strip()
    
        if line == "---------":
            if current_block_idx:
                track_clusters.append(current_block_idx)
                track_times.append(current_block_times)
                current_block_idx = []
                current_block_times = []
            continue

        if 'pass' in line:
            continue

        if "t:" in line:
            cluster_times.append(float(line[2:]))

        try:
            tup = line.split(",")
            trk_idx = int(tup[0])
            trk_time = float(tup[1])
            current_block_idx.append(trk_idx)
            current_block_times.append(trk_time)
        except ValueError:
            # Skip lines like "No valid time" or headers
            continue

    # Catch last block if not followed by dashed line
    if current_block_idx:
        track_clusters.append(current_block_idx)

    print("Parsed time clusters:")
    for cluster in track_clusters:
        print(cluster)
    print("Root script executed successfully.")
except subprocess.CalledProcessError as e:
    print(f"Error executing root script: {e}")
    print(f"Stderr: {e.stderr}")

# flattened_track_indices = list(chain.from_iterable(track_clusters))
# flattened_track_times = [my_branches.Track_time[event_num][t] for i in track_clusters for t in i]

# Collect information for track time histogram
hist_times = []; time_errors = []
pt_wghts = []#; dR_wghts = []; jR_wghts = []; tR_wghts = []; sg_wghts = []; ps_wghts = []
hist_zs = []; z_errors = []
max_pt = -1; max_pt_idx = -1; max_pt_time = -1
hard_scatter_times = []; hard_scatter_zs = []
cluster_zs = []
cluster_z_uncertainties = []
for (num1, cluster) in enumerate(track_clusters):
    this_timecluster = []; this_timeerror = []
    this_ptw = []#; this_dRw = []; this_jRw = []; this_tRw = []; this_sgw = []; this_psw = [];
    this_zcluster = []; this_zerror = []

    zbar_num = 0
    zbar_den = 0
    for (num2,track) in enumerate(cluster):
        this_pt = my_branches.Track_pt[event_num][track]
        this_z0 = my_branches.Track_z0[event_num][track]
        this_var_z0 = my_branches.Track_var_z0[event_num][track]
        this_phi = my_branches.Track_phi[event_num][track]
        this_eta = my_branches.Track_eta[event_num][track]
        this_hasValidTime = my_branches.Track_hasValidTime[event_num][track]

        this_time = track_times[num1][num2]
        this_timeres = my_branches.Track_timeRes[event_num][track]

        zbar_num += this_z0/(this_var_z0)
        zbar_den += 1/(this_var_z0)
        
        if this_pt > max_pt:
            max_pt = this_pt
            max_pt_idx = track
            max_pt_time = this_time

        if (my_branches.Track_truthVtx_idx[event_num][track] == 0):
            hard_scatter_times.append(this_time)
            hard_scatter_zs.append(this_z0)
            
                
        # print(vtx_z + z0)
        this_timecluster.append(this_time); this_timeerror.append(this_timeres)
        this_ptw.append(this_pt)
        this_zcluster.append(this_z0); this_zerror.append(np.sqrt(this_var_z0))
    cluster_zs.append(zbar_num/zbar_den)
    cluster_z_uncertainties.append(1/np.sqrt(zbar_den))
    hist_times.append(this_timecluster); time_errors.append(this_timeerror)
    hist_zs.append(this_zcluster); z_errors.append(this_zerror)
    pt_wghts.append(this_ptw);

weight_cases = {
    # 'Default counts': None,  # No weights, just counts
    'Track pT': pt_wghts,
    '2': pt_wghts,
    # 'exp(-min(dR(track,jet)))': dR_wghts,
    # 'jet pT*exp(-min(dR(track,jet)))': jR_wghts,
    # 'trk pT*exp(-min(dR(track,jet)))': tR_wghts,
    # 'exp(-|s|)': sg_wghts,
    # 'trk pT*exp(-|s|)': ps_wghts,
}

print('--------------------')
print(f'Cluster times: {cluster_times}')
print('--------------------')

print('--------------------')
print(f'Cluster zs: {cluster_zs}')
print(f'Cluster z uncs: {cluster_z_uncertainties}')

print('Cluster dzs:')
for (z,s_z,zs) in zip(cluster_zs,cluster_z_uncertainties,hist_zs):
    regdz = np.abs(reco_hs_z-np.mean(zs))
    print(f'Delta z (mm): {regdz:.2f}')
    
    dz = np.abs(reco_hs_z-z)
    dzres = dz/s_z
    print(f'Delta z (mm, waverage): {dz:.2f}')
    print(f'Delta z (res, waverage): {dzres:.2f}')
print('--------------------')


print('--------------------')
print(f'Hard Scatter Track Times: {hard_scatter_times}')
print('--------------------')

print(f"MAX PT IDX: {max_pt_idx}, PT: {max_pt}, Z0: {my_branches.Track_z0[event_num][max_pt_idx]}")
# find nearest vertex to maxpt track in z
# nearest = 0
# dzNearest = np.abs(my_branches.Track_z0[event_num][max_pt_idx]-my_branches.RecoVtx_z[event_num][nearest])
# for (ivtx, z) in enumerate(my_branches.RecoVtx_z[event_num]):
#     dz = np.abs(my_branches.Track_z0[event_num][max_pt_idx]-z)
#     if dz < dzNearest:
#         nearest=ivtx
#         dzNearest = dz
# nearestz = my_branches.RecoVtx_z[event_num][nearest]

# print(f"NEAREST: {nearest}, nsigma: {dzNearest/np.sqrt(my_branches.Track_var_z0[event_num][max_pt_idx]):.2f}")
# print(f"HARD SCATTER: {0}, nsigma: {np.abs(my_branches.Track_z0[event_num][max_pt_idx]-my_branches.RecoVtx_z[event_num][0])/np.sqrt(my_branches.Track_var_z0[event_num][max_pt_idx]):.2f}")
    
################### Draw Time Histogram #######################
min_time, max_time = min([reco_hs_t, truth_hs_t] + list(chain.from_iterable(hist_times))), max([reco_hs_t, truth_hs_t] + list(chain.from_iterable(hist_times)))
trange = max_time - min_time;
extended_min_time = min_time - 0.05 * trange;
extended_max_time = max_time + 0.05 * trange;

min_z, max_z = min([reco_hs_z, truth_hs_z] + list(chain.from_iterable(hist_zs))), max([reco_hs_z, truth_hs_z] + list(chain.from_iterable(hist_zs)))
zrange = max_z - min_z;
extended_min_z = min_z - 0.05 * zrange;
extended_max_z = max_z + 0.05 * zrange;

filename = f'./figs/nohs/event_display_{file_num}_{event_num:04d}_{vtxID}.pdf'

colors = [
    "#e41a1c",  # red
    "#377eb8",  # blue
    "#4daf4a",  # green
    "#984ea3",  # purple
    "#ff7f00",  # orange
    "#000075",  # navy
    "#bfef45",  # lime
    "#ffff33",  # yellow
    "#aaffc3",  # mint
    "#a65628",  # brown
    "#469990",  # teal
    "#808000",  # olive
    "#dcbeff",  # lavender
    "#f781bf",  # pink
    "#ffd8b1",  # apricot
    "#fffac8",  # beige
    "#999999",  # gray
]


cluster_colors = colors[:len(track_clusters)]
print(cluster_colors)

flattened_hist_times = list(chain.from_iterable(hist_times))
flattened_hist_zs = list(chain.from_iterable(hist_zs))

with PdfPages(filename) as pdf:
    for label, weights in weight_cases.items():
        flattened_weights = list(chain.from_iterable(weights))
        figs, axs = plt.subplots(2,1,figsize=(12, 12))
        time_histo = axs[1]
        event_display = axs[0]

        if '2' in label:
             # Plot histogram with or without weights
            histvals, bin_edges, patches = time_histo.hist(
                hist_times, bins=50, stacked=True, range=(extended_min_time, extended_max_time),
                weights=weights, label='Track Time', lw=1.5, alpha=1.0, color=cluster_colors
            )

            highlight_patch = None
            n_bins = len(bin_edges) - 1
        
            # We'll look for the bin containing max_pt_time
            for bin_idx, (x_left, x_right) in enumerate(zip(bin_edges[:-1], bin_edges[1:])):
                if x_left <= max_pt_time < x_right:
                    if hasattr(histvals, '__len__') and len(histvals) > 0 and hasattr(histvals[0], '__len__'):
                        for container in patches:  # all clusters
                            patch = container.patches[bin_idx]
                            patch.set_hatch('///')
                            patch.set_linewidth(1.5)
                            # Create the legend patch with matching color
                            highlight_patch = mpatches.Patch(
                                facecolor='none',
                                edgecolor='black',
                                hatch='///',
                                linewidth=1.5,
                                label='Max $p_T$ Bin'
                            )
                            break
                    else:
                        patch = patches[bin_idx]
                        patch.set_hatch('///')
                        patch.set_linewidth(1.5)
                        highlight_patch = mpatches.Patch(
                            facecolor='none',
                            edgecolor='black',
                            hatch='///',
                            linewidth=1.5,
                            label='Max $p_T$ Bin'
                        )
                        break
            
            maxval = -1
            if hasattr(histvals, '__len__') and len(histvals) > 0 and hasattr(histvals[0], '__len__'):
                maxval = max(chain.from_iterable(histvals))
            else:
                maxval = histvals.max()

            time_histo.set_ylim(0, 1.1 * maxval)
            truth_text_y = 0.1*maxval
            reco_text_y = 0.1*maxval
            ylabel = 'Track pT'
            cluster_ypoints = [0.1*maxval]*len(cluster_times)
            hist_midpoint = (extended_max_time + extended_min_time) / 2.0
            # Get bin indices where hard scatter times fall
            hs_bin_indices = np.digitize(hard_scatter_times, bin_edges) - 1
            hs_bin_set = set(hs_bin_indices)
        
            # Compute bin centers
            hs_bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

            # Optional: add a marker or text at each bin center
            for i in hs_bin_set:
                center = hs_bin_centers[i]
                time_histo.text(center, 0, '/', ha='center', va='top', fontsize=20, color='blue')
        else:
            for (times, errtime, zs, errz, color) in zip(hist_times, time_errors, hist_zs, z_errors, cluster_colors):
                histvals = time_histo.errorbar(
                    times, zs, fmt='.',
                    xerr=errtime, yerr=errz,
                    lw=1.5, alpha=1.0, color=color,
                    zorder=3
                )
            
            truth_text_y = truth_hs_z
            reco_text_y = reco_hs_z
            cluster_ypoints = cluster_zs
            ylabel = 'Z (mm)'
            time_histo.scatter(hard_scatter_times, hard_scatter_zs,
                               s=20, color='blue', marker='o',label='Hard Scatter Tracks', zorder=10)
            time_histo.plot([time_histo.get_xlim()[0],time_histo.get_xlim()[1]], [reco_hs_z,reco_hs_z],
                            label='Reco Vertex Z', color='black', lw=2, linestyle='--')            
            
        min_y_step = (time_histo.get_ylim()[1] - time_histo.get_ylim()[0])/25
        min_x_step = (time_histo.get_xlim()[1] - time_histo.get_xlim()[0])/50

        time_histo.scatter(cluster_times, cluster_ypoints, marker="*", edgecolors='black',
                           s=500, c=cluster_colors, label='Cluster Positions', zorder=2)
            
        time_histo.annotate("Truth HS Time",
                            xytext=(truth_hs_t-2.5*min_x_step, time_histo.get_ylim()[0]-2*min_y_step),
                            xy=(truth_hs_t, truth_text_y), xycoords='data',
                            ha='left', va='top', color='blue',
                            arrowprops=dict(arrowstyle="->", color='blue',lw=2))

        if (my_branches.RecoVtx_hasValidTime[event_num][0] == 1):
            time_histo.annotate("Reco HS Time",
                                xytext=(reco_hs_t-2.5*min_x_step, time_histo.get_ylim()[0]-3*min_y_step),
                                xy=(reco_hs_t, reco_text_y), xycoords='data',
                                ha='left', va='top', color='black',
                                arrowprops=dict(arrowstyle="->",color='black',lw=2))
        
        if args.extra_time:
            time_histo.annotate("My Time",
                            xytext=(args.extra_time-2.5*min_x_step, time_histo.get_ylim()[0]-4*min_y_step),
                            xy=(args.extra_time, reco_text_y),
                            ha='left', va='top', 
                            color='black',
                            arrowprops=dict(arrowstyle="->", color='black',lw=2))
        
        # # Add any additional lines or markers
        time_histo.axvspan(time_histo.get_xlim()[0], guess_time-3*90, color='grey', alpha=0.2, zorder=1,
                           label=f'Calo Time Exclusion\n[{guess_time-3*90:.2f},{guess_time+3*90:.2f}]')
        time_histo.axvspan(time_histo.get_xlim()[1], guess_time+3*90, color='grey', alpha=0.2, zorder=1)

        # Grab the existing handles and labels
        handles, labels = time_histo.get_legend_handles_labels()

        # Calculate weight sums per cluster
        weight_sums = [len(cluster) for cluster in hist_times] if weights is None else [sum(w) for w in weights]

        # Create dummy legend entries for cluster weight totals
        weight_legend_patches = []
        for i, (color, total_weight) in enumerate(zip(cluster_colors, weight_sums)):
            this_label = f"Σpt = {total_weight:.2f}\ndz={np.abs(reco_hs_z-cluster_zs[i]):.2f}\nt={cluster_times[i]:.2f}"
            patch = mpatches.Patch(facecolor=color, label=this_label)
            weight_legend_patches.append(patch)

        # Add a second legend
        second_legend = time_histo.legend(
            handles=weight_legend_patches,
            loc='upper right',
            title="Cluster Weights",
            fontsize='small',
            title_fontsize='small',
        )
        
        # Add the second legend manually so both show up
        time_histo.add_artist(second_legend)
        time_histo.legend(handles=handles, labels=labels, loc='upper left')
        
        time_histo.set_xlim(extended_min_time, extended_max_time)
        time_histo.set_xlabel('Time (ps)')
        time_histo.set_ylabel(ylabel)
        time_histo.set_title('Z-T Event Display')
        
        ################### Draw Tracks ###############################

        for track in track_info:
            Z = track[0]  # vertex z
            z0 = track[1] # track z0
            x = track[2]  # x length
            y = track[3]  # y length
            status = track[4]
            linestyle = '-'
            if status == 1:
                color = 'blue'  
            elif status == 0:
                color = 'red'
            elif status == 2:
                color = 'blue'
                linestyle = '--'
            elif status == 3:        
                color = 'red'
                linestyle = '--'
            elif status == 4:
                color = 'blue'
                linestyle = '-.'
            elif status == 5:
                color = 'red'
                linestyle = '-.'
            else:
                color = "black"
                # print(z0)
            event_display.plot([Z + z0, Z + z0 + x], [0, y], color=color, linestyle=linestyle)
    
        ################## Draw Reference eta lines ##################################
        def draw_eta_reference_lines(vtx_z, eta_ref, line_length=50, linestyle='dashed'):
            theta_ref_pos = 2 * np.arctan(np.exp(-eta_ref))  # For +eta_ref
            theta_ref_neg = 2 * np.arctan(np.exp(eta_ref))   # For -eta_ref (neg. eta)
            
            # Compute x and y components based on theta
            x_ref_pos = line_length * np.cos(theta_ref_pos)
            y_ref_pos = line_length * np.sin(theta_ref_pos)
            
            x_ref_neg = line_length * np.cos(theta_ref_neg)
            y_ref_neg = -line_length * np.sin(theta_ref_neg)  # Negative y for -eta_ref
            
            # Compute mirrored lines
            y_ref_pos_mirror = -y_ref_pos
            y_ref_neg_mirror = -y_ref_neg

            # Draw the reference lines
            # event_display.plot([vtx_z, vtx_z + x_ref_pos], [0, y_ref_pos], linestyle=linestyle, color='lightgrey',label=fr'$\eta = {eta_ref}$')
            event_display.plot([vtx_z, vtx_z + x_ref_pos], [0, y_ref_pos], linestyle=linestyle, color='lightgrey')
            event_display.plot([vtx_z, vtx_z + x_ref_neg], [0, y_ref_neg], linestyle=linestyle, color='lightgrey')
            event_display.plot([vtx_z, vtx_z + x_ref_pos], [0, y_ref_pos_mirror], linestyle=linestyle, color='lightgrey')
            event_display.plot([vtx_z, vtx_z + x_ref_neg], [0, y_ref_neg_mirror], linestyle=linestyle, color='lightgrey')

        linestyles = ['dotted']
        etas = [2.4]
        for line, eta in zip(linestyles,etas):
            draw_eta_reference_lines(vtx_z, eta_ref=eta, linestyle=line)


        ################## Draw Jets ##################################
        text_index = 0
        # jet_pt        0
        # jet_eta       1
        # jet_phi       2
        # jet_isHS_size 3
        # jet_x         4
        # jet_y         5
        # Rpt           6

        # print(len(jet_info))

        for i in range(len(jet_info)):
            x_jet = jet_info[i][4]  # jet_x is at index 4
            y_jet = jet_info[i][5]  # jet_y is at index 5
            isHS = jet_info[i][3]
            Rpt = jet_info[i][6]
    
            # print("jet#", i, "jet_pt:", jet_info[i][0], "isHS:", isHS, "Rpt:", Rpt)
            
    
            # Condition to filter jets
            color = 'grey'
            if isHS == 1:
                color = 'green'
            elif isHS == 2:
                color = 'blue'
                alpha = 1.0
    
            cone_length = np.sqrt(x_jet**2 + y_jet**2)  # Length of the cone
            cone_width = cone_length * 0.3  # Width proportional to length

            # Perpendicular vector for cone base width
            norm = np.sqrt(x_jet**2 + y_jet**2)
            perp_x = -y_jet / norm * cone_width / 2
            perp_y = x_jet / norm * cone_width / 2

            # Define cone vertices (triangle)
            tip_x, tip_y = vtx_z, 0  # Tip at the vertex z location
            base_left_x, base_left_y = vtx_z + x_jet + perp_x, y_jet + perp_y
            base_right_x, base_right_y = vtx_z + x_jet - perp_x, y_jet - perp_y

            event_display.fill([tip_x, base_left_x, base_right_x], [tip_y, base_left_y, base_right_y], color=color, alpha=0.5)
            jet_text = f"Jet {text_index+1}: $p_T$={jet_info[i][0]:.0f} GeV, $\eta$={jet_info[i][1] :.1f}"
            color2 = 'green' if isHS >= 1 else 'black'  # Set color based on isHS value
            #event_display.text(vtx_z-4.5-0.3, 0.9 - (1.2 + i * 0.1), jet_text, weight='bold', fontsize=12, color=color2)
            event_display.text(vtx_z - 6.5 - 0.3, 0.9 - (1.2 + text_index * 0.1), jet_text,
                               weight='bold', fontsize=12, color=color2)
            text_index += 1
            ################################################################
    

        # Determine the x-coordinate based on the Z coordinate range
        x_coord = vtx_z-6.8  # Set it to the minimum x-coordinate of the line_positions
        y_coord = 0.9  # Set the y-coordinate for the text annotations
        event_display.text(x_coord, y_coord, f"Reco HS (z,t) = ({reco_hs_z:.1f} mm, {reco_hs_t:.1f} ps)", weight='bold', fontsize=12)
        event_display.text(x_coord, y_coord - 0.1, f"Truth HS (z,t) = ({truth_hs_z:.1f} mm, {truth_hs_t:.1f} ps)", weight='bold', fontsize=12)
        # event_display.text(x_coord, y_coord - 0.2, f"Truth Vtx {min_idx} (z,t) = ({matched_z:.1f} mm, {matched_time:.1f} ps)", weight='bold', fontsize=12)
    
        ################## Add tracks and jets legend ###################
        
        legend_entries = [mlines.Line2D([], [], color='blue'),
                          mlines.Line2D([], [], color='red'),
                          mpatches.Rectangle((0, 0), 1, 1,  color='green', alpha=0.5),
                          # mpatches.Rectangle((0, 0), 1, 1,  color='blue', alpha=0.5),
                          mpatches.Rectangle((0, 0), 1, 1,  color='grey', alpha=0.5)]

        label_names = ['Hard Scatter',
                       'Pile-Up',
                       'HS jet',
                       'PU jet']

        eta_legend = [mlines.Line2D([], [], color='lightgrey',linestyle=style) for style,eta in zip(linestyles,etas)]
        eta_labels = [fr'$\eta = {eta}$' for eta in etas]

        # Combine both legends into one
        event_display.legend(legend_entries, label_names, loc='upper left', title='Track and Jet Types', bbox_to_anchor=(0.0, 0.9))

        # Add the legend manually to the plot
        # figs.gca().add_artist(hs_legend)

        ##################### Draw Reco and truth vertices ###############

        # Plot the line for the recovertex_z values
        reco_vertices_z = [z for z, _ in reco_vertices]
        marker_colors_reco = ['blue' if z == vtx_z else 'black' for z in reco_vertices_z]
        event_display.scatter(reco_vertices_z, [0] * len(reco_vertices_z), color=marker_colors_reco, marker='o', s=100)
        event_display.axhline(y=0.0, color='black', linestyle='--')

        truth_vertices_z = [z for z, _ in truth_vertices]
        marker_colors = ['blue' if z == truth_hs_z else
                         'black' for z in truth_vertices_z]
        event_display.scatter(truth_vertices_z, [-0.75] * len(truth_vertices_z), color=marker_colors, marker='|', s=100)
        event_display.axhline(y=-0.75, color='black', linestyle='--')

        label_x = vtx_z + 4.5
        label_y = 0.05
        event_display.text(label_x, label_y, 'Reco vertices', fontsize=12)
        event_display.text(label_x, label_y-0.75, 'Truth vertices', fontsize=12)
        #event_display.text(label_x-3, 0.85, 'ATLAS Simulation Preliminary', fontsize=16, weight='bold', style='italic')
        event_display.text(label_x-2.5, 0.85, 'ATLAS Simulation Internal', fontsize=16, weight='bold', style='italic')


        ################################################################

        event_display.set_ylim(-1.0, 1.0)
        event_display.set_xlim(vtx_z-7.0, vtx_z+7.0)
        event_display.set_title(f'Event# {event_num} : Reco Vertex# {vtxID}')
        event_display.set_xlabel('Z [mm]')
        #event_display.ylabel('R [mm]')
        event_display.set_yticks([])
        # event_display.legend()
        #event_display.grid(True)
        #event_display.savefig(f'figures/fig_{event_num}_{vtxID}_sumpt2.png')
        # filename = 
                
        # Save current figure page
        pdf.savefig(figs)
        plt.close(figs)  # Close figure to free memory
        print(f"Event display page printed to {filename}")


