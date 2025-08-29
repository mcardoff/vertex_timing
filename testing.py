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


events_to_try = [86]
event_num = 404
file_num = '000007'
root_file = uproot.open(f'../ntuple-hgtd/user.mcardiff.45809429.Output._{file_num}.SuperNtuple.root') 
tree = root_file["ntuple"]
branch = tree.arrays([
    'TruthVtx_z',
    'TruthVtx_time',
    'TruthVtx_isHS',
    'RecoVtx_isHS',
    'RecoVtx_x',
    'RecoVtx_y',
    'RecoVtx_z',
    'RecoVtx_time',
    'RecoVtx_hasValidTime',
    'Track_truthVtx_idx',
    'RecoVtx_sumPt2',
    'RecoVtx_track_idx',
    'AntiKt4EMTopoJets_ghostTrack_idx',
    'Track_z0',
    'Track_var_z0',
    'Track_d0',
    'Track_var_d0',
    'Track_pt',
    'Track_eta',
    'Track_time',
    'Track_timeRes',
    'Track_hasValidTime',
    'Track_qOverP',
    'Track_theta',
    'Track_phi',
    'Track_charge',
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

track_clusters = []
cluster_times = []
track_times = []
try:
    macro_call = f'runHGTD_Clustering.cxx("{file_num}",{event_num},false,true)'
    result = subprocess.run(['root', '-l', '-q', '-b', macro_call],
                            check=True, capture_output=True, text=True)
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
        except ValueError:
            # Skip lines like "No valid time" or headers
            continue

    # Catch last block if not followed by dashed line
    if current_block_idx:
        track_clusters.append(current_block_idx)

    print("Parsed time clusters:")
    # for (cluster,times) in zip(track_clusters,track_times):
    #     print(cluster)
    #     print(times)
    print("Root script executed successfully.")
except subprocess.CalledProcessError as e:
    print(f"Error executing root script: {e}")
    print(f"Stderr: {e.stderr}")

waverage_num = 0
waverage_den = 0
    
for (i, cluster) in enumerate(track_clusters):
    for (j,track) in enumerate(cluster):
        if (branch.Track_truthVtx_idx[event_num][track] == 0 and i != 0):
            particle = branch.Track_truthPart_idx[event_num][track]
            delta = -1
            if particle != -1:
                delta = np.abs(track_times[i][j] - branch.TruthPart_prodVtx_time[event_num][particle])
            print(f'Track t: {track_times[i][j]:.2f}+/-{branch.Track_timeRes[event_num][track]:.2f}')
            print(f'Track pT: {branch.Track_pt[event_num][track]:.2f}')
            print(f'Delta t truth {delta:.2f}')
            print('--------------------')
            waverage_num += np.random.normal(branch.RecoVtx_time[event_num][0],30)/(branch.Track_timeRes[event_num][track])**2
            waverage_den += 1/((branch.Track_timeRes[event_num][track])**2)

# truth_hs_vtx_z = branch.TruthVtx_z[event_num][0]
# reco_hs_vtx_z = branch.RecoVtx_z[event_num][0]
# print(f'Reco z: {reco_hs_vtx_z}')
# print(f'Truth z: {truth_hs_vtx_z}')
# close_to_truthhs_idx = 0
# min_delta = 10000000000
# for (i,z) in enumerate(branch.RecoVtx_z[event_num]):
#     delta = np.abs(z-truth_hs_vtx_z)
#     if (delta < min_delta):
#         close_to_truthhs_idx = i
#         min_delta = delta

# print(f'Closest Reco Vtx, #{close_to_truthhs_idx} z: {branch.RecoVtx_z[event_num][close_to_truthhs_idx]}')

for (t,track) in zip(track_times[0],track_clusters[0]):
    print(f'Track t: {t:.2f}+/-{branch.Track_timeRes[event_num][track]:.2f}')
    print('--------------------')
    waverage_num += np.random.normal(branch.RecoVtx_time[event_num][0],30)/(branch.Track_timeRes[event_num][track])**2
    waverage_den += 1/((branch.Track_timeRes[event_num][track])**2)


print(f'weighted_avg: {waverage_num/waverage_den}')
