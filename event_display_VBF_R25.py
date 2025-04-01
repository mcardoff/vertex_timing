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
import awkward as ak
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

#############################################################
import argparse
parser = argparse.ArgumentParser(description='Process event and vertex data from ROOT file.')
parser.add_argument('--event_num', type=int, required=True, help='Event number to process')
parser.add_argument('--vtxID', type=int, required=True, help='Vertex ID to process')
parser.add_argument('--file_num', type=str, required=True, help='Specific ID of Filen to process')
args = parser.parse_args()

event_num = args.event_num
vtxID = args.vtxID
file_num = args.file_num
#############################################################

#root_file = uproot.open('CaloTimingNtuple_ttbarSingleLep_FPT_n100_test_20241211.root')
#root_file = uproot.open('user.scheong.42871997.Output._000050.SuperNtuple.root')
root_file = uproot.open(f'./ntuple/user.scheong.42871997.Output._{file_num}.SuperNtuple.root')
tree = root_file["ntuple"]

my_branches = tree.arrays([
    'TruthVtx_z',
    'TruthVtx_time',
    'TruthVtx_isHS',
    'RecoVtx_isHS',
    'RecoVtx_z',
    'RecoVtx_time',
    'Track_truthVtx_idx',
    'RecoVtx_sumPt2',
    # 'RecoVtx_isPU',
    'RecoVtx_track_idx',
    'AntiKt4EMTopoJets_track_idx',
    'Track_z0',
    'Track_pt',
    'Track_time',
    'Track_hasValidTime',
    'Track_qOverP',
    'Track_theta',
    'Track_phi',
    'Track_var_z0',
    'AntiKt4EMTopoJets_pt',
    'AntiKt4EMTopoJets_eta',
    'AntiKt4EMTopoJets_phi',
    'AntiKt4EMTopoJets_truthHSJet_idx',
    'TruthHSJet_pt',
    'TruthHSJet_eta',
    'TruthHSJet_phi',
])

# my_branches = tree.arrays(['TruthVtx_z', 'RecoVtx_sumPt2', 'RecoVtx_z'])

vtx_z = my_branches.RecoVtx_z[event_num][vtxID]
vtx_t = my_branches.RecoVtx_time[event_num][vtxID]
sumpt = my_branches.RecoVtx_sumPt2[event_num][vtxID]
truth_z = my_branches.TruthVtx_z[event_num][0]
truth_t = my_branches.TruthVtx_time[event_num][0]

selected_HS_vtx_id = None
max_sumpt = float('-inf')
my_track_z0=[]


# Iterate through the vertex z-values and sumpt values
for idx, z in enumerate(my_branches.RecoVtx_z[event_num]):
    if abs(z - truth_z) <= 5:
        sumpt_value = my_branches.RecoVtx_sumPt2[event_num][idx]
        #print("reco vertex within 5 mm of truth z : ", z, "vtx id:", idx, "sumpt :", sumpt_value, "|(z - truth_z)| :", abs(z - truth_z))
        if  sumpt_value > max_sumpt:
            selected_HS_vtx_id = idx
            max_sumpt = my_branches.RecoVtx_sumPt2[event_num][idx]
            
print(f"python event_display_VBF_R25.py --event_num {event_num} --vtxID {selected_HS_vtx_id}")


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
print("isHS =",my_branches.RecoVtx_isHS[event_num][vtxID]
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
    if abs(z - vtx_z) <= 5:
        truth_vertices.append((z, idx))
        
#####################################################

# Get all tracks within 3 sigma of selected vertex
# connected_tracks = my_branches.RecoVtx_track_idx[event_num][vtxID]
connected_tracks = []
for idx in range(len(my_branches.Track_z0[event_num])):
    dz = my_branches.Track_z0[event_num][idx] - my_branches.RecoVtx_z[event_num][vtxID]
    nsigma = dz / sqrt(my_branches.Track_var_z0[event_num][idx])
    isHS = my_branches.Track_truthVtx_idx[event_num][idx] == 0
    if(abs(nsigma) < 3.0 and my_branches.Track_pt[event_num][idx] > 1.0 and isHS):
        #   and my_branches.Track_hasValidTime[event_num][idx] == 1
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
    for j in range(len(my_branches.AntiKt4EMTopoJets_track_idx[event_num])):
        if my_branches.AntiKt4EMTopoJets_pt[event_num][j] < 30.0:
            continue
        # njets_b1=njets_b1+1    
        trackPT0 = 0
        trackPT = 0
        jet_pt_05_rpt = 0
        jet_pt_01_rpt = 0

        for q in range(len(my_branches.AntiKt4EMTopoJets_track_idx[event_num][j])):
            # Track selection cuts
            jdx = my_branches.AntiKt4EMTopoJets_track_idx[event_num][j][q]
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
        
    
    
    # new_sumpt= new_sumpt + (track_pT ** 2)

    pz = track_pT * math.sinh(track_eta)
    signX = track_eta / abs(track_eta)
    signY = math.sin(track_phi) / abs(math.sin(track_phi))
    theta = math.atan(track_pT / abs(pz))
    x = (track_pT / 2) * math.cos(theta) * signX
    y = (track_pT / 2) * math.sin(theta) * signY

    Track_truthVtx_id = my_branches.Track_truthVtx_idx[event_num][idx]
    if (Track_truthVtx_id!=0):
        continue

    status = my_branches.Track_hasValidTime[event_num][idx]
    my_branches.TruthVtx_isHS[event_num][Track_truthVtx_id]
    
    if (my_branches.TruthVtx_isHS[event_num][Track_truthVtx_id]==1):
        num_HS_tracks=num_HS_tracks+1
        track_info.append([vtx_z, z0, x, y, status])
    

    
print("number of HS tracks : ", num_HS_tracks)

for truth_jet_idx in range(len(my_branches.TruthHSJet_pt[event_num])):
    truth_jet_pt  = my_branches.TruthHSJet_pt[event_num][truth_jet_idx]
    if truth_jet_pt < 30.0:
        continue
    truth_jet_eta = my_branches.TruthHSJet_eta[event_num][truth_jet_idx]
    truth_jet_phi = my_branches.TruthHSJet_phi[event_num][truth_jet_idx]
    truth_jet_pz = truth_jet_pt * np.sinh(truth_jet_eta)

    truth_jet_signX = np.sign(truth_jet_eta) if truth_jet_eta != 0 else 1  # Avoid division by zero
    truth_jet_signY = np.sign(np.sin(truth_jet_phi)) if np.sin(truth_jet_phi) != 0 else 1

    truth_jet_theta = np.arctan(truth_jet_pt / abs(truth_jet_pz))
    truth_jet_x = (truth_jet_pt / 40) * np.cos(truth_jet_theta) * truth_jet_signX
    truth_jet_y = (truth_jet_pt / 40) * np.sin(truth_jet_theta) * truth_jet_signY
    truth_jet_tuple = (
        truth_jet_pt, truth_jet_eta, truth_jet_phi, 2, truth_jet_x, truth_jet_y, 1.0)

    jet_info.append(truth_jet_tuple)


# Collect information for track time histogram
hist_times = []
for (track, z0) in enumerate(my_branches.Track_z0[event_num]):
    this_pt = my_branches.Track_pt[event_num][track]
    this_hasValidTime = my_branches.Track_hasValidTime[event_num][track]
    if (this_hasValidTime == 0) or this_pt < 1.0:
        continue

    nsigma = (z0 - vtx_z)/np.sqrt(my_branches.Track_var_z0[event_num][track])
    if np.abs(nsigma) > 3.0:
        continue

    hist_times.append(my_branches.Track_time[event_num][track])

# Plot the graph
figs, axs = plt.subplots(2,1,figsize=(12, 12))
time_histo = axs[1]
event_display = axs[0]

################### Draw Time Histogram #######################
histvals = time_histo.hist(hist_times,50,label='Track Time', lw=5, alpha=0.35, color='blue')
time_histo.plot([vtx_t,vtx_t],[0, histvals[0].max()], 'green', label='Reco Vtx Time')
time_histo.plot([truth_t,truth_t],[0, histvals[0].max()], 'red', label='Truth Vtx Time')
time_histo.legend()
time_histo.set_ylim(0, 1+histvals[0].max())
time_histo.set_xlabel("Track Time [ps]")
time_histo.set_ylabel("Frequency")
time_histo.set_title("Track Time Histogram")

################### Draw Tracks ###############################

for track in track_info:
    Z = track[0]  # vertex z
    z0 = track[1] # track z0
    x = track[2]  # x length
    y = track[3]  # y length
    status = track[4]
    if status == 1:
        color = 'blue'  
    elif status == 0:
        color = 'red'
    elif status == 2:        
        color = 'green'
    elif status == 3:        
        color = 'cyan' 
    else:
        color = 'black' 
    event_display.plot([Z + z0, Z + z0 + x], [0, y], color=color)
    
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
    
    # Define label position
    # label_x_eta = vtx_z + 3.5
    # label_y_eta = 0.05

    # Adjust y position for multiple labels to avoid overlap
    # offset = 0.5 * (eta_ref - 2.9)  # Slight shift based on eta_ref value
    # event_display.plot([label_x_eta - 0.4, label_x_eta], [label_y_eta - 0.90 - offset, label_y_eta - 0.90 - offset], linestyle=linestyle, color='lightgrey', linewidth=1.5)
    # event_display.text(label_x_eta + 0.2, label_y_eta - 0.90 - offset, , fontsize=12, color='lightgrey')

linestyles = ['dashed', 'dotted']
etas = [2.4, 2.0]
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
    
    print("jet#", i, "jet_pt:", jet_info[i][0], "isHS:", isHS, "Rpt:", Rpt)

    
    # Condition to filter jets
    if num_HS_tracks > 0 and abs(vtx_z - truth_z) < 2:
        if isHS == 0:
            continue  # Skip non-HS jets when this condition is met
    else:
        if isHS == 0 and Rpt < 0.02:
            continue

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
    event_display.text(vtx_z - 4.5 - 0.3, 0.9 - (1.2 + text_index * 0.1), jet_text,
             weight='bold', fontsize=12, color=color2)
    text_index += 1
################################################################
    
# Determine the x-coordinate based on the Z coordinate range
x_coord = vtx_z-4.8  # Set it to the minimum x-coordinate of the line_positions
y_coord = 0.9  # Set the y-coordinate for the text annotations
event_display.text(x_coord, y_coord, f"Reco (z,t) = ({vtx_z:.1f} mm, {vtx_t:.1f} ps)", weight='bold', fontsize=12)
event_display.text(x_coord, y_coord - 0.1, f"Truth (z,t) = ({truth_z:.1f} mm, {truth_t:.1f} ps)", weight='bold', fontsize=12)
    
################## Add tracks and jets legend ###################

legend_entries = [mlines.Line2D([], [], color='blue'),
                  mlines.Line2D([], [], color='red'),
                  mpatches.Rectangle((0, 0), 1, 1,  color='green', alpha=0.5),
                  mpatches.Rectangle((0, 0), 1, 1,  color='blue', alpha=0.5),
                  mpatches.Rectangle((0, 0), 1, 1,  color='grey', alpha=0.5)]

label_names = ['Valid Time',
               'No Valid Time',
               'HS jet',
               'Truth jet',
               'PU jet']

eta_legend = [mlines.Line2D([], [], color='lightgrey',linestyle=style) for style,eta in zip(linestyles,etas)]
eta_labels = [fr'$\eta = {eta}$' for eta in etas]

# Combine both legends into one
event_display.legend(legend_entries+eta_legend, label_names+eta_labels, loc='upper right', title='Track and Jet Types', bbox_to_anchor=(1.0, 0.9))

# Add the legend manually to the plot
# figs.gca().add_artist(hs_legend)

##################### Draw Reco and truth vertices ###############

# Plot the line for the recovertex_z values
reco_vertices_z = [z for z, _ in reco_vertices]
marker_colors_reco = ['blue' if z == vtx_z else 'black' for z in reco_vertices_z]
event_display.scatter(reco_vertices_z, [0] * len(reco_vertices_z), color=marker_colors_reco, marker='o', s=100)
event_display.axhline(y=0.0, color='black', linestyle='--')

truth_vertices_z = [z for z, _ in truth_vertices]
marker_colors = ['blue' if z == truth_z else 'black' for z in truth_vertices_z]
event_display.scatter(truth_vertices_z, [-0.75] * len(truth_vertices_z), color=marker_colors, marker='|', s=100)
event_display.axhline(y=-0.75, color='black', linestyle='--')

label_x = vtx_z + 3.5
label_y = 0.05
event_display.text(label_x, label_y, 'Reco vertices', fontsize=12)
event_display.text(label_x, label_y-0.75, 'Truth vertices', fontsize=12)
#event_display.text(label_x-3, 0.85, 'ATLAS Simulation Preliminary', fontsize=16, weight='bold', style='italic')
event_display.text(label_x-2.2, 0.85, 'ATLAS Simulation Internal', fontsize=16, weight='bold', style='italic')


################################################################

event_display.set_ylim(-1.0, 1.0)
event_display.set_xlim(vtx_z-5.0, vtx_z+5.0)
event_display.set_title(f'Event# {event_num} : Vertex# {vtxID}')
event_display.set_xlabel('Z [mm]')
#event_display.ylabel('R [mm]')
event_display.set_yticks([])
# event_display.legend()
#event_display.grid(True)
#event_display.savefig(f'figures/fig_{event_num}_{vtxID}_sumpt2.png')
plt.savefig(f'./figs/trackhists/{file_num}/tjeventdisplay_{event_num}_{vtxID}.pdf')

print(f"Event display has been saved as ./figs/trackhists/{file_num}/tjeventdisplay_{event_num}_{vtxID}.pdf")

#plt.show()



