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

python -m Martinoid --seq "Na-Na"