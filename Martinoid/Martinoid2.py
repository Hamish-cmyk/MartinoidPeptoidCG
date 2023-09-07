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


    
class Martinoid:
    def build(self):
        for seq_i, residue in enumerate(self.seq):
            #######################################
            # Add Beads from Data.csv
            last_i = self.beads.index.shape[0]
            new_beads = ParseData(residue, "beads")
            new_index = np.array(new_beads.index) + last_i
            new_beads = new_beads.set_index(new_index)
            new_beads["residue"] = [self.residue_n]*new_beads.shape[0]
            new_beads["resname"] = [residue]*new_beads.shape[0]
            self.beads = pandas.concat((self.beads, new_beads))
            self.res_length.append(len(new_beads))
            #######################################
            # Add BB-BB bond --> Data.csv "BB = []"
            previous_BB_residue = max([seq_i - 1, 0])
            previous_BB_residue = self.seq[previous_BB_residue]
            if np.where(self.beads["i"]=="BB")[0].shape[0] >= 2:  
                bond = np.where(self.beads["i"]=="BB")[0][-2:]
                self.bonds.loc[self.bonds.index.shape[0]] = [bond[0], bond[1], ParseData(previous_BB_residue, "BB")[0], ParseData(previous_BB_residue, "BB")[1]]
            #######################################
            # Add BB-BB-BB Angle --> Data.csv "AA = []"
            # "AL" is angle linear | "AH" is angle helical
            if np.where(self.beads["i"]=="BB")[0].shape[0] >= 3:
                BB = np.where(self.beads["i"]=="BB")[0][-3:]
                if self.Helical == True:
                    self.angles.loc[self.angles.index.shape[0]] = [BB[0], BB[1], BB[2], ParseData(self.seq[seq_i-2], "AH")[0],ParseData(previous_BB_residue, "AH")[1], "BB-BB-BB"]
                else:
                    self.angles.loc[self.angles.index.shape[0]] = [BB[0], BB[1], BB[2], ParseData(previous_BB_residue, "AL")[0],ParseData(previous_BB_residue, "AL")[1], "BB-BB-BB"]
            #######################################
            # Add BB-BB-BB-BB dihedral
            # "DL" is dihedral linear | "DH" is dihedral helical
            if np.where(self.beads["i"]=="BB")[0].shape[0] >= 4:
                if self.Helical == True:
                    BB = np.where(self.beads["i"]=="BB")[0][-4:]
                    self.dihedrals.loc[self.dihedrals.index.shape[0]] = [BB[0], BB[1], BB[2], BB[3], ParseData(self.seq[seq_i-3], "DH")[0],ParseData(self.seq[seq_i-3], "DH")[1], "BB-BB-BB-BB"]
                else:
                    BB = np.where(self.beads["i"] == "BB")[0][-4:]
                    self.dihedrals.loc[self.dihedrals.index.shape[0]] = [BB[0], BB[1], BB[2], BB[3], ParseData(previous_BB_residue, "DL")[0],ParseData(previous_BB_residue, "DL")[1], "BB-BB-BB-BB"]
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
                    self.constraints.loc[self.constraints.index.shape[0]] = new_bonds.loc[bond][list("ijL")]
                else:
                    self.bonds.loc[self.bonds.index.shape[0]] = new_bonds.loc[bond]
            ############################################################################################################
            new_angles = ParseData(residue, "angles")
            new_angles["i"] += last_i
            new_angles["j"] += last_i
            new_angles["k"] += last_i
            new_angles["comment"] = "BB-SC1-SC2"
            for angle in new_angles.index:
                self.angles.loc[self.angles.index.shape[0]] = new_angles.loc[angle]
            ############################################################################################################
            new_dihedrals       = ParseData(residue, "dihedrals")                   ## add sum to all columns for that species
            new_dihedrals["i"] += last_i
            new_dihedrals["j"] += last_i
            new_dihedrals["k"] += last_i
            new_dihedrals["l"] += last_i
            new_dihedrals["comment"] = "BB-SC2-SC3-SC1"
            for dihed in new_dihedrals.index:                                       ## loop em in
                self.dihedrals.loc[self.dihedrals.index.shape[0]] = new_dihedrals.loc[dihed]
            ############################################################################################################
            self.residue_n += 1
            
    def set_angles(self):
        #####################################################################################################
        ## NOTE: now setting SC-BB-BB angles after dataframe build loop - this way we can capture all of them
        ## and allows discrimination of the final residue angle which is set to be acute for reasons disucssed
        ## in the model manuscript.
        #####################################################################################################
        backbone = []
        if self.residue_n > 1:
            for ind in self.beads.index:
                beadID = self.beads['i'].iloc[ind]
                if beadID == 'BB':
                    backbone.append(ind)
            # pos1, pos2, pos3 = 0,1,2 unused variables?
            for each in range(0, len(backbone)):
                centre = backbone[each]
                if self.res_length[each] > 1:
                    if each != (len(backbone)-1):
                        self.angles.loc[self.angles.index.shape[0]] = [centre+1,centre,backbone[each+1],'110','20',"SC1-BB-BB"]
                    else:
                        self.angles.loc[self.angles.index.shape[0]] = [centre+1,centre,backbone[each-1],'65','80',"SC1-BB-BB"]
                else:
                    pass
    
    def build_positions(self):
        ###############################################################################
        #### Generate CG bead positions --> mdtraj loop through 'beads' DataFrame #####
        ###############################################################################
        self.bonds["i"] = self.bonds["i"].astype(np.uint64)
        self.bonds["j"] = self.bonds["j"].astype(np.uint64)
        ###############################################################################
        ## Single Mononmer Case
        if len(self.seq) == 1:
            if self.NTC:
                self.beads.at[self.beads[self.beads["i"] == "BB"].index[0], "type"] = "Qd"
            else:
                self.beads.at[self.beads[self.beads["i"] == "BB"].index[0], "type"] = "Nda"
        ## Dimers and Beyond
        else:
            ####################################################################
            if self.NTC:
                self.beads.at[self.beads[self.beads["i"] == "BB"].index[0], "type"] = "Qd"
            else:
                if self.NTCC:
                    self.beads.at[self.beads[self.beads["i"] == "BB"].index[0], "type"] = "Na"     ## Acylated NTER
                else:
                    self.beads.at[self.beads[self.beads["i"] == "BB"].index[0], "type"] = "Nda"    ## Neutral NTER
            ####################################################################
            if self.CTC:
                self.beads.at[self.beads[self.beads["i"] == "BB"].index[-1], "type"] = "Nda"       ## Amidated CTER
            else:
                if self.CTCC:
                    self.beads.at[self.beads[self.beads["i"] == "BB"].index[-1], "type"] = "Qa"    ## Charged CTER
                else:
                    self.beads.at[self.beads[self.beads["i"] == "BB"].index[-1], "type"] = "P3"    ## Neutral CTER
            #####################################################################
        template = mdtraj.core.topology.Topology()
        coords = np.array([[0.0, 0.0, 0.0]])
        coords_beads = ["BB"]
        index = 0
        sign = 1
        while index <= self.beads["i"][1:].shape[0]:
            bead = self.beads.at[index, "i"]
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
                SC_i = int(self.beads.at[index, "i"][-1])-1
                template.add_atom(f"SC{SC_i}", mdtraj.core.element.carbon, res)
                new_side_chain = SideChainCoords[self.beads.at[index, "resname"]] + new_pos
                new_side_chain[:,1]*=sign
                coords = np.vstack((coords, new_side_chain[SC_i]))
                coords_beads += ["SC"] 
            index += 1
        traj = mdtraj.Trajectory(coords, template)
        traj.save_pdb(f"{self.name}.pdb")
        print(f"Writing Structure >>> {self.name}.pdb")
        ### Set the variables as internal class parameters for debugging
        self.coords = coords
        
    def post_process_pdb(self):
        pdbcontent = readin(f"{self.name}.pdb")
        with open(f"{self.name}.pdb", 'w') as opdb:
            for line in pdbcontent.split("\n"):
                if "ATOM" in line:
                    opdb.write(line)
                    opdb.write("\n")
    
    def debug_positions(self):
        ##############################
        #### SIMPLE DEBUG of beads ###
        ##############################
        local_i = 0
        for i in range(self.beads.shape[0]):
            if self.beads.at[i, "i"] == "BB":
                c = "blue"
            else:
                c = "red"
            if self.beads.at[i, "i"] == "BB":
                local_i = 0
            plt.scatter([coords[i,0]], [coords[i,1]], color=c)
            plt.text(coords[i,0], coords[i,1], str(local_i))
            local_i += 1
        
    def write_itp(self):
        ###########################
        #### ITP WRITER IS HERE ###
        ###########################
        itp = open(f"{self.name}.itp", 'w')
        command = " ".join(sys.argv)
        itp.write(f'; input >> {command}\n')
        itp.write(f"""[ moleculetype ]
; Name         Exclusions
Peptoid   1\n
""")
        itp.write("[ atoms ]\n")
        for bead in self.beads.index:
            i = bead+1
            typ = self.beads.at[bead, "type"]
            B = self.beads.at[bead, "i"]
            resname = self.beads.at[bead, "resname"]
            ###################
            if resname == 'Na':
                resname = 'Nx'
            ###################
            resnum = self.beads.at[bead, "residue"] + 1
            charge = charges[typ]
            writ = f"\t{i}\t{typ}\t{resnum}\t{resname}\t{B}\t {i}  {charge} ; \n"
            print(writ)
            itp.write(f"\t{i}\t{typ}\t{resnum}\t{resname}\t{B}\t {i}  {charge} ; \n")
        itp.write("\n[ bonds ]\n")
        for bond in self.bonds.index:
            i = int(self.bonds.at[bond, "i"]+1)
            j = int(self.bonds.at[bond, "j"]+1)
            L = self.bonds.at[bond, "L"]
            k = self.bonds.at[bond, "k"]
            b1 = self.beads.at[i-1, "i"]
            b2 = self.beads.at[j-1, "i"]
            itp.write(f"\t{i}\t{j}\t1\t{L}\t\t{k} ; {b1}-{b2}\n")
        itp.write("\n[ constraints ]\n")
        for constraint in self.constraints.index:
            i = int(self.constraints.at[constraint, "i"]+1)
            j = int(self.constraints.at[constraint, "j"]+1)
            L = self.constraints.at[constraint, "L"]
            b1 = self.beads.at[i-1, "i"]
            b2 = self.beads.at[j-1, "i"]
            itp.write(f"\t{i}\t{j}\t1\t{L} ; {b1}-{b2}\n")
        itp.write("\n[ angles ]\n")
        exclude = []
        for angle in self.angles.index:
            i = int(self.angles.at[angle, "i"]+1)
            j = int(self.angles.at[angle, "j"]+1)
            k = int(self.angles.at[angle, "k"]+1)
            Angle = self.angles.at[angle, "angle"]
            FC = self.angles.at[angle, "FC"]
            #############################
            type = '2'
            comment = self.angles.at[angle, "comment"]
            if Angle == 'X' and FC == 'Y':
                Angle  = '0'
                FC     = '1'
                type   = '8'
                ignore = f'{i} {k}'
                exclude.append(ignore)
                comment = 'DW_Potential'
            #############################
            itp.write(f"\t{i}\t{j}\t{k}\t{type}\t{Angle}\t{FC} ; {comment}\n")
        if self.dihedrals.shape[0] > 0:
            itp.write("\n[ dihedrals ]\n")
            for dihe in self.dihedrals.index:
                i = int(self.dihedrals.at[dihe, "i"]+1)
                j = int(self.dihedrals.at[dihe, "j"]+1)
                k = int(self.dihedrals.at[dihe, "k"]+1)
                l = int(self.dihedrals.at[dihe, "l"]+1)
                dihedral = self.dihedrals.at[dihe, "dihedral"]
                FC = self.dihedrals.at[dihe, "FC"]
                if self.dihedrals.at[dihe, "comment"] == "BB-BB-BB-BB":
                    if self.Helical == True:
                        comment = self.dihedrals.at[dihe, "comment"]
                        print(f"\t{i}\t{j}\t{k}\t{l}\t1\t{dihedral}\t{FC}\t1 ; {comment}")
                        itp.write(f"\t{i}\t{j}\t{k}\t{l}\t1\t{dihedral}\t{FC}\t1 ; {comment} \n")
                    else:
                        comment = self.dihedrals.at[dihe, "comment"]
                        itp.write(f"\t{i}\t{j}\t{k}\t{l}\t1\t{dihedral}\t{FC}\t1 ; {comment} \n")
                else:
                    comment = self.dihedrals.at[dihe, "comment"]
                    itp.write(f"\t{i}\t{j}\t{k}\t{l}\t2\t{dihedral}\t{FC} ; {comment} \n")
        itp.write("\n[ exclusions ]\n")
        for ex in range(0,len(exclude)):
            itp.write(f"\t{exclude[ex]}\n")
        itp.close()
        print(f"Writing itp file  >>> {self.name}.itp")
        print("\nMartinoid for to say:", np.random.choice(TPB_quotes),'\n')
                
    def __init__(self, sequence, NTC=False, CTC=False, NTCC=False, CTCC=False, Linear=False, Helical=False):
        # Unpack arguments
        self.seq = sequence.split("-")
        self.name = name = "".join([f"{x.replace('0', '')}-" for x in sequence.split("-")])[:-1]
        self.NTC = NTC
        self.CTC = CTC
        self.NTCC = NTCC
        self.CTCC = CTCC
        self.Linear = Linear
        self.Helical = Helical
        print("Name:", name)
        # Setup dataframes
        ##############################################################################################
        self.beads       = pandas.DataFrame(columns=["i", "type", "residue", "resname"])
        self.bonds       = pandas.DataFrame(columns=["i", "j", "L", "k"], dtype=np.float64)
        self.constraints = pandas.DataFrame(columns=["i", "j", "L"], dtype=np.float64)
        self.angles      = pandas.DataFrame(columns=["i", "j", "k", "angle", "FC", "comment"])
        self.dihedrals   = pandas.DataFrame(columns=["i", "j", "k", "l", "dihedral", "FC", "comment"])
        ##############################################################################################
        self.residue_n  = 0
        self.res_length = []
        # Build the peptoid topology
        self.build()
        self.set_angles()
        # guess some reasonable positions for the beads
        self.build_positions()
        self.post_process_pdb() # clean up the pdb file
        #self.debug_positions()
        # Write the itp topology
        self.write_itp()

