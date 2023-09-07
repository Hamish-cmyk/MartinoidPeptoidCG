# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 18:22:22 2021

@author: avtei
"""
from Martinoid._martinoid_data import *
import numpy as np
import pandas, sys, time, argparse
import warnings
import mdtraj
from matplotlib import pyplot as plt, colors
warnings.filterwarnings("ignore")
np.random.seed(int(time.time()))

DEBUG   = True
Linear  = False
Helical = False



if len(sys.argv) < 3:
    if '-h' in sys.argv or "--help" in sys.argv:

        print(statement)
        sys.exit()
    else:
        print("Input Data Required...")
        sys.exit()
else:
    resid     = sys.argv[1]
    NTC  = True
    CTC  = True
    NTCC = False
    CTCC = False
    if len(sys.argv) > 2:
        #########################
        if "-nt" in sys.argv:
            NTC  = False
        if "-ntc" in sys.argv:
            NTC  = False
            NTCC = True
        #########################
        if "-ct" in sys.argv:
            CTC  = False
        if "-ctc" in sys.argv:
            CTC  = False
            CTCC = True
        #########################
        if "Linear" in sys.argv:        ## apply loose 180 BB-BB-BB-BB dihedral
            Linear = True
        if "Helical" in sys.argv:       ## apply rigid 90 degree BB-BB-BB-BB dihedral
            Helical = True
        #########################
seq = resid.split("-")
name = "".join([f"{x.replace('0', '')}-" for x in resid.split("-")])[:-1]
print("\nName:", name)
print(NTC)
##############################################################################################
beads       = pandas.DataFrame(columns=["i", "type", "residue", "resname"])
bonds       = pandas.DataFrame(columns=["i", "j", "L", "k"], dtype=np.float64)
constraints = pandas.DataFrame(columns=["i", "j", "L"], dtype=np.float64)
angles      = pandas.DataFrame(columns=["i", "j", "k", "angle", "FC", "comment"])
dihedrals   = pandas.DataFrame(columns=["i", "j", "k", "l", "dihedral", "FC", "comment"])
##############################################################################################
residue_n  = 0
res_length = []
for seq_i, residue in enumerate(seq):
    #######################################
    # Add Beads from Data.csv
    last_i = beads.index.shape[0]
    new_beads = ParseData(residue, "beads")
    new_index = np.array(new_beads.index) + last_i
    new_beads = new_beads.set_index(new_index)
    new_beads["residue"] = [residue_n]*new_beads.shape[0]
    new_beads["resname"] = [residue]*new_beads.shape[0]
    beads = pandas.concat((beads, new_beads))
    res_length.append(len(new_beads))
    #######################################
    # Add BB-BB bond --> Data.csv "BB = []"
    previous_BB_residue = max([seq_i - 1, 0])
    previous_BB_residue = seq[previous_BB_residue]
    if np.where(beads["i"]=="BB")[0].shape[0] >= 2:  
        bond = np.where(beads["i"]=="BB")[0][-2:]
        bonds.loc[bonds.index.shape[0]] = [bond[0], bond[1], ParseData(previous_BB_residue, "BB")[0], ParseData(previous_BB_residue, "BB")[1]]
    #######################################
    # Add BB-BB-BB Angle --> Data.csv "AA = []"
    # "AL" is angle linear | "AH" is angle helical
    if np.where(beads["i"]=="BB")[0].shape[0] >= 3:
        BB = np.where(beads["i"]=="BB")[0][-3:]
        if Helical == True:
            angles.loc[angles.index.shape[0]] = [BB[0], BB[1], BB[2], ParseData(seq[seq_i-2], "AH")[0],ParseData(previous_BB_residue, "AH")[1], "BB-BB-BB"]
        else:
            angles.loc[angles.index.shape[0]] = [BB[0], BB[1], BB[2], ParseData(previous_BB_residue, "AL")[0],ParseData(previous_BB_residue, "AL")[1], "BB-BB-BB"]
    #######################################
    # Add BB-BB-BB-BB dihedral
    # "DL" is dihedral linear | "DH" is dihedral helical
    if np.where(beads["i"]=="BB")[0].shape[0] >= 4:
        if Helical == True:
            BB = np.where(beads["i"]=="BB")[0][-4:]
            dihedrals.loc[dihedrals.index.shape[0]] = [BB[0], BB[1], BB[2], BB[3], ParseData(seq[seq_i-3], "DH")[0],ParseData(seq[seq_i-3], "DH")[1], "BB-BB-BB-BB"]
        else:
            BB = np.where(beads["i"] == "BB")[0][-4:]
            dihedrals.loc[dihedrals.index.shape[0]] = [BB[0], BB[1], BB[2], BB[3], ParseData(previous_BB_residue, "DL")[0],ParseData(previous_BB_residue, "DL")[1], "BB-BB-BB-BB"]
            pass
    #############################################################################################################
    # MAKE SIDECHAIN COMPONENTS
    ############################################################################################################
    new_bonds = ParseData(residue, "bonds")
    new_bonds["i"] += last_i
    new_bonds["j"] += last_i
    for bond in new_bonds.index:
        isConstraint = str(new_bonds.loc[bond]["k"])
        if "constraint" in isConstraint.lower():
            constraints.loc[constraints.index.shape[0]] = new_bonds.loc[bond][list("ijL")]
        else:
            bonds.loc[bonds.index.shape[0]] = new_bonds.loc[bond]
    ############################################################################################################
    new_angles = ParseData(residue, "angles")
    new_angles["i"] += last_i
    new_angles["j"] += last_i
    new_angles["k"] += last_i
    new_angles["comment"] = "BB-SC1-SC2"
    for angle in new_angles.index:
        angles.loc[angles.index.shape[0]] = new_angles.loc[angle]
    ############################################################################################################
    new_dihedrals       = ParseData(residue, "dihedrals")                   ## add sum to all columns for that species
    new_dihedrals["i"] += last_i
    new_dihedrals["j"] += last_i
    new_dihedrals["k"] += last_i
    new_dihedrals["l"] += last_i
    new_dihedrals["comment"] = "BB-SC2-SC3-SC1"
    for dihed in new_dihedrals.index:                                       ## loop em in
        dihedrals.loc[dihedrals.index.shape[0]] = new_dihedrals.loc[dihed]
    ############################################################################################################
    residue_n += 1
#####################################################################################################
## NOTE: now setting SC-BB-BB angles after dataframe build loop - this way we can capture all of them
## and allows discrimination of the final residue angle which is set to be acute for reasons disucssed
## in the model manuscript.
#####################################################################################################
backbone = []
if residue_n > 1:
    for ind in beads.index:
        beadID = beads['i'].iloc[ind]
        if beadID == 'BB':
            backbone.append(ind)
    pos1, pos2, pos3 = 0,1,2
    for each in range(0, len(backbone)):
        centre = backbone[each]
        if res_length[each] > 1:
            if each != (len(backbone)-1):
                angles.loc[angles.index.shape[0]] = [centre+1,centre,backbone[each+1],'110','20',"SC1-BB-BB"]
            else:
                angles.loc[angles.index.shape[0]] = [centre+1,centre,backbone[each-1],'65','80',"SC1-BB-BB"]
        else:
            pass
###############################################################################
#### Generate CG bead positions --> mdtraj loop through 'beads' DataFrame #####
###############################################################################
bonds["i"] = bonds["i"].astype(np.uint64)
bonds["j"] = bonds["j"].astype(np.uint64)
###############################################################################
## Single Mononmer Case
if len(seq) == 1:
    if NTC:
        beads.at[beads[beads["i"] == "BB"].index[0], "type"] = "Qd"
    else:
        beads.at[beads[beads["i"] == "BB"].index[0], "type"] = "Nda"
## Dimers and Beyond
else:
    print(NTC,NTCC)
    ####################################################################
    if NTC:
        beads.at[beads[beads["i"] == "BB"].index[0], "type"] = "Qd"
    else:
        if NTCC:
            beads.at[beads[beads["i"] == "BB"].index[0], "type"] = "Na"     ## Acylated NTER
        else:
            beads.at[beads[beads["i"] == "BB"].index[0], "type"] = "Nda"    ## Neutral NTER
    ####################################################################
    if CTC:
        beads.at[beads[beads["i"] == "BB"].index[-1], "type"] = "Nda"       ## Amidated CTER
    else:
        if CTCC:
            beads.at[beads[beads["i"] == "BB"].index[-1], "type"] = "Qa"    ## Charged CTER
        else:
            beads.at[beads[beads["i"] == "BB"].index[-1], "type"] = "P3"    ## Neutral CTER
    #####################################################################
template = mdtraj.core.topology.Topology()
coords = np.array([[0.0, 0.0, 0.0]])
coords_beads = ["BB"]
index = 0
sign = 1
while index <= beads["i"][1:].shape[0]:
    bead = beads.at[index, "i"]
    new_pos = coords[np.where(np.array(coords_beads) == "BB")[0][-1]].copy()
    new_pos[2] = np.random.random()/10
    if bead == "BB":
        current_chain = template.add_chain()
        res = template.add_residue("W1000", current_chain)
        template.add_atom("BB", mdtraj.core.element.oxygen, res)
        new_pos[0] += 0.4
        new_pos[1] = 0
        if index > 0:
            coords = np.vstack((coords, new_pos))
            coords_beads.append("BB")
        if sign == 1:
            sign = -1
        else:
            sign = 1
    else:
        SC_i = int(beads.at[index, "i"][-1])-1
        template.add_atom(f"SC{SC_i}", mdtraj.core.element.carbon, res)
        new_side_chain = SideChainCoords[beads.at[index, "resname"]] + new_pos
        new_side_chain[:,1]*=sign
        coords = np.vstack((coords, new_side_chain[SC_i]))
        coords_beads += ["SC"] 
    index += 1
traj = mdtraj.Trajectory(coords, template)
traj.save_pdb(f"{name}.pdb")
print(f"Writing Structure >>> {name}.pdb")
##############################
#### SIMPLE DEBUG of beads ###
##############################
local_i = 0
if DEBUG:
    for i in range(beads.shape[0]):
        if beads.at[i, "i"] == "BB":
            c = "blue"
        else:
            c = "red"
        if beads.at[i, "i"] == "BB":
            local_i = 0
        plt.scatter([coords[i,0]], [coords[i,1]], color=c)
        plt.text(coords[i,0], coords[i,1], str(local_i))
        local_i += 1
###########################
#### ITP WRITER IS HERE ###
###########################
itp = open(f"{name}.itp", 'w')
NAME = f"{name}"
command = " ".join(sys.argv)
itp.write(f'; input >> {command}\n')
itp.write(f"""[ moleculetype ]
; Name         Exclusions
Peptoid   1\n
""")
itp.write("[ atoms ]\n")
for bead in beads.index:
    i = bead+1
    typ = beads.at[bead, "type"]
    B = beads.at[bead, "i"]
    resname = beads.at[bead, "resname"]
    ###################
    if resname == 'Na':
        resname = 'Nx'
    ###################
    resnum = beads.at[bead, "residue"] + 1
    charge = charges[typ]
    writ = f"\t{i}\t{typ}\t{resnum}\t{resname}\t{B}\t {i}  {charge} ; \n"
    print(writ)
    itp.write(f"\t{i}\t{typ}\t{resnum}\t{resname}\t{B}\t {i}  {charge} ; \n")
itp.write("\n[ bonds ]\n")
for bond in bonds.index:
    i = int(bonds.at[bond, "i"]+1)
    j = int(bonds.at[bond, "j"]+1)
    L = bonds.at[bond, "L"]
    k = bonds.at[bond, "k"]
    b1 = beads.at[i-1, "i"]
    b2 = beads.at[j-1, "i"]
    itp.write(f"\t{i}\t{j}\t1\t{L}\t\t{k} ; {b1}-{b2}\n")
itp.write("\n[ constraints ]\n")
for constraint in constraints.index:
    i = int(constraints.at[constraint, "i"]+1)
    j = int(constraints.at[constraint, "j"]+1)
    L = constraints.at[constraint, "L"]
    b1 = beads.at[i-1, "i"]
    b2 = beads.at[j-1, "i"]
    itp.write(f"\t{i}\t{j}\t1\t{L} ; {b1}-{b2}\n")
itp.write("\n[ angles ]\n")
exclude = []
for angle in angles.index:
    i = int(angles.at[angle, "i"]+1)
    j = int(angles.at[angle, "j"]+1)
    k = int(angles.at[angle, "k"]+1)
    Angle = angles.at[angle, "angle"]
    FC = angles.at[angle, "FC"]
    #############################
    type = '2'
    comment = angles.at[angle, "comment"]
    if Angle == 'X' and FC == 'Y':
        Angle  = '0'
        FC     = '1'
        type   = '8'
        ignore = f'{i} {k}'
        exclude.append(ignore)
        comment = 'DW_Potential'
    #############################
    itp.write(f"\t{i}\t{j}\t{k}\t{type}\t{Angle}\t{FC} ; {comment}\n")
if dihedrals.shape[0] > 0:
    itp.write("\n[ dihedrals ]\n")
    for dihe in dihedrals.index:
        i = int(dihedrals.at[dihe, "i"]+1)
        j = int(dihedrals.at[dihe, "j"]+1)
        k = int(dihedrals.at[dihe, "k"]+1)
        l = int(dihedrals.at[dihe, "l"]+1)
        dihedral = dihedrals.at[dihe, "dihedral"]
        FC = dihedrals.at[dihe, "FC"]
        if dihedrals.at[dihe, "comment"] == "BB-BB-BB-BB":
            if Helical == True:
                comment = dihedrals.at[dihe, "comment"]
                print(f"\t{i}\t{j}\t{k}\t{l}\t1\t{dihedral}\t{FC}\t1 ; {comment}")
                itp.write(f"\t{i}\t{j}\t{k}\t{l}\t1\t{dihedral}\t{FC}\t1 ; {comment} \n")
            else:
                comment = dihedrals.at[dihe, "comment"]
                itp.write(f"\t{i}\t{j}\t{k}\t{l}\t1\t{dihedral}\t{FC}\t1 ; {comment} \n")
        else:
            comment = dihedrals.at[dihe, "comment"]
            itp.write(f"\t{i}\t{j}\t{k}\t{l}\t2\t{dihedral}\t{FC} ; {comment} \n")
itp.write("\n[ exclusions ]\n")
for ex in range(0,len(exclude)):
    itp.write(f"\t{exclude[ex]}\n")
itp.close()
print(f"Writing itp file  >>> {name}.itp")
print("\nMartinoid for to say:", np.random.choice(TPB_quotes),'\n')
sys.exit()