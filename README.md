# MartinoidPeptoidCG

This program was originally written as part of the publication title "..", published in the ***Journal of Very Cool Science***

This module was inspired by martinize (http://cgmartini.nl/index.php/tools2/proteins-and-bilayers/204-martinize) and has been created to perform automatic topology building of peptoids within the MARTINI forcefield (v2.1) in the GROMACS program.

A key difference between Martinoid and Martinize is that the former does not require an input all-atom peptoid structure while Martinize does. This has obvious advantages but does mean the output guessed CG structures [generally] require a greater degree of minimization.

## Installation
(Will be availible on PyPi later)

	git clone https://github.com/Hamish-cmyk/MartinoidPeptoidCG
	cd MartinoidPeptoidCG
	pip install .

## Usage

The primary argument passed is the peptoid sequence, this must be in the correct format according to the Glasgow convention (https://doi.org/10.17868/strath.00085559).

Martinoid can either be used directly from the terminal once installed:
For example:

	python -m Martinoid --seq "Na-Nt-Nfe"

Which will output two files: Na-Nt-Nfe.itp & Na-Nt-Nfe.pdb

The topology (Na-Nt-Nfe.itp) looks like:
```
[ moleculetype ]
; Name         Exclusions
Peptoid   1

[ atoms ]
	1	Qd	1	Nx	BB	 1  1.0 ; 
	2	Na	2	Nt	BB	 2  0.0 ; 
	3	SP1	2	Nt	SC1	 3  0.0 ; 
	4	Nda	3	Nfe	BB	 4  0.0 ; 
	5	SC1	3	Nfe	SC1	 5  0.0 ; 
	6	SC5	3	Nfe	SC2	 6  0.0 ; 
	7	SC5	3	Nfe	SC3	 7  0.0 ; 
	8	SC5	3	Nfe	SC4	 8  0.0 ; 

[ bonds ]
	1	2	1	0.34		8515.0 ; BB-BB
	2	3	1	0.333		6627.0 ; BB-SC1
	2	4	1	0.34		8515.0 ; BB-BB
	4	5	1	0.298		7638 ; BB-SC1
	5	6	1	0.234		9500 ; SC1-SC2

[ constraints ]
	6	7	1	0.27 ; SC2-SC3
	7	8	1	0.27 ; SC3-SC4
	8	6	1	0.27 ; SC4-SC2

[ angles ]
	1	2	4	2	138.0	50.0 ; BB-BB-BB
	4	5	6	2	147	37 ; BB-SC1-SC2
	5	6	7	8	0	1 ; DW_Potential
	5	6	8	8	0	1 ; DW_Potential
	3	2	4	2	110	20 ; SC1-BB-BB
	5	4	2	2	65	80 ; SC1-BB-BB

[ dihedrals ]
	5	7	8	6	2	0	50 ; BB-SC2-SC3-SC1 

[ exclusions ]
	5 7
	5 8
}
```

and the coordinates (Na-Nt-Nfe.pdb) has a guess structure which looks like:
![alt text](images/Na-Nt-Nfe.png "Na-Nt-Nfe.pdb")
