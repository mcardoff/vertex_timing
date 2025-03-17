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
    'TruthVtx_isHS',
    'RecoVtx_isHS',
    'Track_truthVtx_idx',
    'TruthVtx_isHS',
    'RecoVtx_sumPt2',
    # 'RecoVtx_isPU',
    'RecoVtx_track_idx',
    'AntiKt4EMTopoJets_track_idx',
    'RecoVtx_z',
    'Track_z0',
    'Track_pt',
    'Track_hasValidTime',
    'Track_qOverP',
    'Track_theta',
    'Track_phi',
    'Track_var_z0',
    'AntiKt4EMTopoJets_pt',
    'AntiKt4EMTopoJets_eta',
    'AntiKt4EMTopoJets_phi',
    'AntiKt4EMTopoJets_truthHSJet_idx'
])

# my_branches = tree.arrays(['TruthVtx_z', 'RecoVtx_sumPt2', 'RecoVtx_z'])

vtx_z = my_branches.RecoVtx_z[event_num][vtxID]
sumpt = my_branches.RecoVtx_sumPt2[event_num][vtxID]
truth_z = my_branches.TruthVtx_z[event_num][0]

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
    if(abs(nsigma) < 2.0 and my_branches.Track_pt[event_num][idx] > 1.0 and my_branches.Track_hasValidTime[event_num][idx] == 1):
        connected_tracks.append(idx)

track_info = []
jet_info = []
num_jets = 0
new_sumpt=0
sum_pt_W =0

num_HS_tracks=0

for idx in connected_tracks:
    track_z0 = my_branches.Track_z0[event_num][idx]
    p = abs(1 / (my_branches.Track_qOverP[event_num][idx]))
    track_eta = -np.log(math.tan((my_branches.Track_theta[event_num][idx]) / 2))
    track_pT = (p / (np.cosh(track_eta))) / 1000
    track_phi = my_branches.Track_phi[event_num][idx]
    z0 = track_z0 - vtx_z
    
    ################################
    #if (track_pT>25):
    #    continue
    ################################
    
    my_track_z0.append(track_z0)

    ####################################################
    ### Calculate minDR ################################
    ####################################################
    
    minDr = 1000;
    closestJetIndex = -1;
    jpt_val = 0.0;

    num_jet_match_01=0
    num_jet_match_05=0
    
    njets_b1 = 0
    njets_b2 = 0

    # Loop over jets
    for j in range(len(my_branches.AntiKt4EMTopoJets_track_idx[event_num])):
        if my_branches.AntiKt4EMTopoJets_pt[event_num][j] < 30.0:
            continue
        njets_b1=njets_b1+1
    
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
                continue
            trackPT += pt2
    
        Rpt = trackPT / my_branches.AntiKt4EMTopoJets_pt[event_num][j]
        #print("jet#", j, "jet_pt", my_branches.AntiKt4EMTopoJets_pt[event_num][j], "Rpt:", Rpt)
            
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
        
        #num_jets += 1  # Increment counter
        
        if Rpt < 0.02:
        #if Rpt < 0.00:
            continue
            
    ####################################################
    ####################################################
        njets_b2=njets_b2+1

        deta = my_branches.AntiKt4EMTopoJets_eta[event_num][j] - track_eta
        dphi = my_branches.AntiKt4EMTopoJets_phi[event_num][j] - my_branches.Track_phi[event_num][idx]
        
        if dphi > pi:
            dphi -= 2 * pi
        Dr = sqrt(deta**2 + dphi**2)
    
        # Check if the current jet is closer than the previously found closest jet
        if Dr < minDr and Dr <= 0.8:
        #if Dr < minDr:
            minDr = Dr  # Update minimum Dr
            closestJetIndex = j  # Update index of closest jet
            jpt = my_branches.AntiKt4EMTopoJets_pt[event_num][closestJetIndex]
            jpt_val = jpt
    
        jet_pt_01_rpt += my_branches.AntiKt4EMTopoJets_pt[event_num][j]
        num_jet_match_01 += 1
    
        if Rpt > 0.5:
            jet_pt_05_rpt += my_branches.AntiKt4EMTopoJets_pt[event_num][j]
            num_jet_match_05 += 1
            
    #print(njets_b1, njets_b2)
    
    if minDr == 1000:
        minDr = 1.0

    #print("minDr :", minDr)
    #print("jpt :", jpt_val)

    pt_W = (track_pT ** 2) * (jpt_val ** 2) / minDr
    sum_pt_W = sum_pt_W+pt_W

    #print("pt_W :", pt_W)
    
    ####################################################
    ####################################################
    
    
    new_sumpt= new_sumpt + (track_pT ** 2)

    pz = track_pT * math.sinh(track_eta)
    signX = track_eta / abs(track_eta)
    signY = math.sin(track_phi) / abs(math.sin(track_phi))
    theta = math.atan(track_pT / abs(pz))
    x = (track_pT / 2) * math.cos(theta) * signX
    y = (track_pT / 2) * math.sin(theta) * signY

    #status = my_branches.track_status[event_num][idx]
    #status = my_branches.Track_isTruthHS[event_num][idx]
    Track_truthVtx_id = my_branches.Track_truthVtx_idx[event_num][idx]
    if (Track_truthVtx_id==-1):
        continue

    status = my_branches.TruthVtx_isHS[event_num][Track_truthVtx_id]
    
    if (status==1):
        num_HS_tracks=num_HS_tracks+1
    
    track_info.append([vtx_z, z0, x, y, status])
    #track_info.append(compute_track_line(track_pT, track_eta, track_phi, vtx_z, z0, status))

    
#print("jet_pt_list size: ", len(jet_pt_list))

print("number of HS tracks : ", num_HS_tracks)

# Plot the graph
plt.figure(figsize=(12, 6))


################### Draw Tracks ###############################

for track in track_info:
    Z = track[0]  # vertex z
    z0 = track[1] # track z0
    x = track[2]  # x length
    y = track[3]  # y length
    status = track[4]
    if status == 1:
        color = 'blue'  # Set color to blue for hard-scatter tracks
    elif status == 0:
        color = 'red'
    elif status == 2:        
        color = 'green'
    elif status == 3:        
        color = 'cyan' 
    else:
        color = 'black' 
    plt.plot([Z + z0, Z + z0 + x], [0, y], color=color)
    
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
    plt.plot([vtx_z, vtx_z + x_ref_pos], [0, y_ref_pos], linestyle=linestyle, color='lightgrey',label=fr'$\eta = {eta_ref}$')
    plt.plot([vtx_z, vtx_z + x_ref_neg], [0, y_ref_neg], linestyle=linestyle, color='lightgrey')
    plt.plot([vtx_z, vtx_z + x_ref_pos], [0, y_ref_pos_mirror], linestyle=linestyle, color='lightgrey')
    plt.plot([vtx_z, vtx_z + x_ref_neg], [0, y_ref_neg_mirror], linestyle=linestyle, color='lightgrey')
    
    # Define label position
    # label_x_eta = vtx_z + 3.5
    # label_y_eta = 0.05

    # Adjust y position for multiple labels to avoid overlap
    # offset = 0.5 * (eta_ref - 2.9)  # Slight shift based on eta_ref value
    # plt.plot([label_x_eta - 0.4, label_x_eta], [label_y_eta - 0.90 - offset, label_y_eta - 0.90 - offset], linestyle=linestyle, color='lightgrey', linewidth=1.5)
    # plt.text(label_x_eta + 0.2, label_y_eta - 0.90 - offset, , fontsize=12, color='lightgrey')

linestyles = ['dashed', 'dotted']
etas = [2.4, 2.0]
for line, eta in zip(linestyles,etas):
    draw_eta_reference_lines(vtx_z, eta_ref=eta, linestyle=line)


################## Draw Jets ##################################
text_index = 0

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

    color = 'green' if isHS >= 1 else 'grey'
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

    plt.fill([tip_x, base_left_x, base_right_x], [tip_y, base_left_y, base_right_y], color=color, alpha=0.5)
    jet_text = f"Jet {text_index+1}: $p_T$={jet_info[i][0]:.0f} GeV, $\eta$={jet_info[i][1] :.1f}"
    color2 = 'green' if isHS >= 1 else 'black'  # Set color based on isHS value
    #plt.text(vtx_z-4.5-0.3, 0.9 - (1.2 + i * 0.1), jet_text, weight='bold', fontsize=12, color=color2)
    plt.text(vtx_z - 4.5 - 0.3, 0.9 - (1.2 + text_index * 0.1), jet_text,
             weight='bold', fontsize=12, color=color2)
    text_index += 1
################################################################
    
# Determine the x-coordinate based on the Z coordinate range
x_coord = vtx_z-4.8  # Set it to the minimum x-coordinate of the line_positions
y_coord = 0.9  # Set the y-coordinate for the text annotations
plt.text(x_coord, y_coord, f"Reco z = {vtx_z:.1f} mm", weight='bold', fontsize=12)
plt.text(x_coord, y_coord - 0.1, f"Truth z = {truth_z:.1f} mm", weight='bold', fontsize=12)
plt.text(x_coord, y_coord - 0.2, f"Sum $\\mathbf{{p_T^2}}$ = {new_sumpt:.2e} $\\mathbf{{GeV^2}}$", weight='bold', fontsize=12)
plt.text(x_coord, y_coord - 0.3, f"Sum $\\mathbf{{p_TW}}$ = {sum_pt_W:.2e} $\\mathbf{{GeV^4}}$", weight='bold', fontsize=12)
    
################## Add tracks and jets legend ###################

track_legend = [mlines.Line2D([], [], color='blue', label='HS tracks'),
                mlines.Line2D([], [], color='red', label='PU tracks')]

jet_legend = [mpatches.Rectangle((0, 0), 1, 1,  color='green', alpha=0.5, label='HS jet'),
              mpatches.Rectangle((0, 0), 1, 1,  color='grey', alpha=0.5, label='PU jet')]

eta_legend = [mlines.Line2D([], [], color='lightgrey',linestyle=style, label=fr'$\eta = {eta}$') for style,eta in zip(linestyles,etas)]

# Combine both legends into one
hs_legend = plt.legend(handles=track_legend + jet_legend, loc='upper right', title='Track and Jet Types', bbox_to_anchor=(1.0, 0.9))

# Add the legend manually to the plot
plt.gca().add_artist(hs_legend)

##################### Draw Reco and truth vertices ###############

# Plot the line for the recovertex_z values
reco_vertices_z = [z for z, _ in reco_vertices]
marker_colors_reco = ['blue' if z == vtx_z else 'black' for z in reco_vertices_z]
plt.scatter(reco_vertices_z, [0] * len(reco_vertices_z), color=marker_colors_reco, marker='o', s=100)
plt.axhline(y=0.0, color='black', linestyle='--')

truth_vertices_z = [z for z, _ in truth_vertices]
marker_colors = ['blue' if z == truth_z else 'black' for z in truth_vertices_z]
plt.scatter(truth_vertices_z, [-0.75] * len(truth_vertices_z), color=marker_colors, marker='|', s=100)
plt.axhline(y=-0.75, color='black', linestyle='--')

label_x = vtx_z + 3.5
label_y = 0.05
plt.text(label_x, label_y, 'Reco vertices', fontsize=12)
plt.text(label_x, label_y-0.75, 'Truth vertices', fontsize=12)
#plt.text(label_x-3, 0.85, 'ATLAS Simulation Preliminary', fontsize=16, weight='bold', style='italic')
plt.text(label_x-2.2, 0.85, 'ATLAS Simulation Internal', fontsize=16, weight='bold', style='italic')


################################################################

plt.ylim(-1.0, 1.0)
plt.xlim(vtx_z-5.0, vtx_z+5.0)
plt.title(f'Event# {event_num} : Vertex# {vtxID}')
plt.xlabel('Z [mm]')
#plt.ylabel('R [mm]')
plt.yticks([])
# plt.legend()
#plt.grid(True)
#plt.savefig(f'figures/fig_{event_num}_{vtxID}_sumpt2.png')
plt.savefig(f'./eventdisplays/0fjet/{file_num}/fig_{event_num}_{vtxID}.pdf')

print(f"Event display has been saved as figures/fig_{event_num}_{vtxID}.pdf")

#plt.show()



