"""
EVENT DISPLAY CODE
CREATED BY: wasikul.islam@cern.ch
MODIFIED FOR USE BY: mcardiff@brandeis.edu
"""

import argparse
import random
import subprocess
from itertools import chain

import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

# --- Argument Parsing and Data Loading ---
# parser = argparse.ArgumentParser(description='Process event and vertex data.')
# parser.add_argument('--event_num', type=int, required=True)
# parser.add_argument('--file_num', type=str, required=True)
# parser.add_argument('--extra_time', type=float, required=False)
# args = parser.parse_args()

event_num, file_num = 599, '000002'

IDEALEFF = False
# filename = f'./VDT_TALK_2/RPT/event_display_{file_num}_{event_num:04d}.pdf'
# colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628", "#f781bf", "#999999"]

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
    'Track_recoVtx_idx',
    'Track_truthVtx_idx',
    # reco jet properties
    'AntiKt4EMTopoJets_pt' ,
    'AntiKt4EMTopoJets_eta',
    'AntiKt4EMTopoJets_phi',
    'AntiKt4EMTopoJets_truthHSJet_idx',
    'TruthHSJet_pt',
    'TruthHSJet_eta',
    'TruthITPUJet_pt',
    'TruthITPUJet_eta',
    'TruthOOTPUJet_pt',
    'TruthOOTPUJet_eta',
    'Track_nearestVtx_sig',
    'Track_nearestVtx_idx',
])

print('='*50)
print('Jet Information: Truth HS Jets')
for (jet_pt, jet_eta) in zip(branch.TruthHSJet_pt[event_num], branch.TruthHSJet_eta[event_num]):
    print(f'Jet pT:  {jet_pt :.2f}')
    print(f'Jet eta: {jet_eta:.2f}')
    print('-'*10)

print('Jet Information: Truth In Time PU Jets')
for (jet_pt, jet_eta) in zip(branch.TruthITPUJet_pt[event_num], branch.TruthITPUJet_eta[event_num]):
    print(f'Jet pT:  {jet_pt :.2f}')
    print(f'Jet eta: {jet_eta:.2f}')
    print('-'*10)

print('Jet Information: Truth Out of Time PU Jets')
for (jet_pt, jet_eta) in zip(branch.TruthOOTPUJet_pt[event_num], branch.TruthOOTPUJet_eta[event_num]):
    print(f'Jet pT:  {jet_pt :.2f}')
    print(f'Jet eta: {jet_eta:.2f}')
    print('-'*10)
print('='*50)

print('Jet Information: AK4EMTopo Reco Jets')
for (jet_pt, jet_eta, truthjets) in zip(branch.AntiKt4EMTopoJets_pt[event_num],
                                        branch.AntiKt4EMTopoJets_eta[event_num],
                                        branch.AntiKt4EMTopoJets_truthHSJet_idx[event_num]):
    print(f'Jet pT:  {jet_pt :.2f}')
    print(f'Jet eta: {jet_eta:.2f}')
    print(f'Associated Truth HS Jets: {truthjets}')
    print('-'*10)
print('='*50)

reco_hs_z = branch.RecoVtx_z[event_num][0]
reco_hs_t = branch.RecoVtx_time[event_num][0]

truth_hs_z = branch.TruthVtx_z[event_num][0]
truth_hs_t = branch.TruthVtx_time[event_num][0]

# find closest associated vertex for each track
connected_tracks = []
for (idx, eta) in enumerate(branch.Track_eta[event_num]):
    in_hgtd = abs(eta) > 2.38 and abs(eta) < 4.0
    time_valid = True if IDEALEFF else branch.Track_hasValidTime[event_num][idx] == 1
    track_quality = branch.Track_quality[event_num][idx] == 1
    dz = branch.Track_z0[event_num][idx] - branch.RecoVtx_z[event_num][0]
    nsigma = abs(dz / np.sqrt(branch.Track_var_z0[event_num][idx]))
    if nsigma < 3.0 and branch.Track_pt[event_num][idx] > 1.0 and in_hgtd and time_valid:
        connected_tracks.append(idx)

print(f'Found {len(connected_tracks)} tracks')

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
            
    if closest_reco_vtx >= 0:
        trk_reco_vtx = branch.Track_recoVtx_idx[event_num][idx_trk]
        trk_truth_vtx = branch.Track_truthVtx_idx[event_num][idx_trk]
        ntuple_nearest_sig = branch.Track_nearestVtx_sig[event_num][idx_trk]
        ntuple_nearest_vtx = int(branch.Track_nearestVtx_idx[event_num][idx_trk])

        trk_nearest_z = branch.RecoVtx_z[event_num][closest_reco_vtx]
        ntuple_nearest_z = branch.RecoVtx_z[event_num][ntuple_nearest_vtx]


    eta_trk = branch.Track_eta[event_num][idx_trk]
    pt_trk = branch.Track_pt[event_num][idx_trk]
        
    print(50*'=')
    print(f'Track {idx_trk}')
    print(f'pt: {pt_trk:.2f}, z0: {z0_trk:.2f}, eta: {eta_trk:.2f}')
    print(f'RecoVtx Association: {trk_reco_vtx}')
    print(f'TruthVtx Association: {trk_truth_vtx}')
    print(f'Calculated Closest Vertex (nsigma = {closest_nsigma:.2f}): {closest_reco_vtx} vtx_z={trk_nearest_z:.2f}')
    print(f'NTuple Closest Vertex (nsigma = {ntuple_nearest_sig:.2f}): {ntuple_nearest_vtx} vtx_z={ntuple_nearest_z:.2f}')
    print(50*'=')
    
