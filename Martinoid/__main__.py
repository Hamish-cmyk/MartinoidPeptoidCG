# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 11:02:13 2023

@author: Alexander van Teijlingen
"""

# Here is the behavior we want it to execute if it is run as a standalone program, ie. like python -m Martinioid ...

import argparse
from Martinoid.Martinoid2 import *

statement = '''
Martinoid2.py is a script to create Coarse Grain Martinoid 
input files of peptoids, ready for use in the molecular 
dynamics simulations package Gromacs.\n
Primary input/output
--------------------
Martinoid2.py [sequence] [Flags] [Secondary Structure]\n\n
Sequence             -- Written using the Glasgow convention (https://doi.org/10.17868/strath.00085559)
                        with a '-' spacer between each molecules (Example: tri sarcosine would be 'Na-Na-Na').\n\n
Flags                 - Martinoid allows users to explicitly set termini backbone charge states:\n
                        -nt  : is a neutral N- terminus with hydrogen bonding acceptor and donor character.
                        -ntc : is an acetylated N-terminus with only hydrogen bonding acceptor character.
                        -ct  : is a neutral carboxylate C-terminus.
                        -ctc : is a charged carboxylate C-terminus with acceptor character\n
                        If no flag is provided then the script defaults to a charged donor N-terminus (Qd) and a
                        non-polar amidated C-terminus (Nda) which reflects many cases in the literature.\n\n
Secondary Structure --  Martinoid version 1.0 has two secondary structure options for sequences longer than trimers: 
                        'Linear' and 'Helical'. The former has a loose linear backbone dihedral while the latter has
                        rigid angle/dihedrals to enforce a Right Handed (RH) helix. Additional secondary structures 
                        are planned for future versions.\n
Output              --  Martinoid will produce an *itp and *pdb file for the given structure. The *itp file is ready 
                        to use with the molecule name "Peptoid" which can be used as a t-coupling group. The *pdb
                        structure will require minimization/equilibration prior to use, though is adequate for system
                        construction.\n 
Directory Structure
----------------   
To function correctly "Martinoid2.py" must be in the directory above ~/Data and in the same directory as "_martinoid_data.py".
We have found it useful to make sub-dirs within a "Martinoid" directory, generate *itp files in these and then reference to these 
from elsewhere in the file system. 
                
Additional Notes
----------------     
The popular chiral residues (Nfes/Nfer), for which handedness is lost through coarse graining, are represented by a single monomer 
'Nfex' which reproduces their non-bonded character through bead type selection. The sarcosine monomer, Na, has the resname 'Nx' because
Gromacs interprets 'Na' as an ion and then proceeds to cause drama.\n

Happy CG Peptoid Simulating!\n\n\t(Hamish W. A. Swanson, hamish.swanson@strath.ac.uk, 05.09.2023)
----------------
'''



if __name__ == "__main__":
    print(statement)
    parser = argparse.ArgumentParser(description="")
    
    parser.add_argument('--seq', required=True, metavar="Sequence", type=str, help='Peptoid sequence, according to the Glasgow convention (https://doi.org/10.17868/strath.00085559)')
    parser.add_argument('-nt', action='store_true', help='Neutral N- terminus with hydrogen bonding acceptor and donor character')
    parser.add_argument('-ntc', action='store_true', help='Acetylated N-terminus with only hydrogen bonding acceptor character')
    parser.add_argument('-ct', action='store_true', help='Neutral carboxylate C-terminus')
    parser.add_argument('-ctc', action='store_true', help='Charged carboxylate C-terminus with acceptor character')
    parser.add_argument('--linear', action='store_true', help='Apply loose 180 BB-BB-BB-BB dihedral')
    parser.add_argument('--helical', action='store_true', help='Apply rigid 90 degree BB-BB-BB-BB dihedral')
    
    
    args = parser.parse_args()
    
    # Check we have the required arguments
    assert isinstance(args.seq, str), "A valid peptoid sequence must be specified with the --seq parameter, such as --seq Na-Na-Na" + statement
    Sequence = args.seq
    
    #Set variables from the optional arguements
    NTC  = True if args.nt == False else False
    CTC  = True if args.ct == False else False
    NTCC = args.ntc
    CTCC = args.ctc
    Linear = args.linear
    Helical = args.helical
    # check that they haven't specified linear and helical together
    assert int(Linear)+int(Helical < 2), "You cannot species --linear and --helical, its one, the other or neither."

    #print(args)
    
    peptoid = Martinoid(Sequence, NTC, CTC, NTCC, CTCC, Linear, Helical)
    
