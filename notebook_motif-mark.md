Motif Mark Project
==================

*********** 2022-02-15 ************

## Setting up pycairo

	$ brew install pkg-config
	$ brew install cairo
	$ pip3 install pycairo

	Successfully installed pycairo-1.20.1

## Run program

	$ python3 motif-mark-oop.py -f Figure_1.fasta -m Fig_1_motifs.txt 

## Jason's suggestions for error handling

raise Exception('Error 505: Silly guy')

## Other notes

- exons are uppercase
- motifs doesn't matter if it's upper or lower case
- motifs can be rev comp
- it's been reverse complemented already (but be sure to say it's verse complimented in the final figure)
- these are from UCSC database
- any functions associated with an object should be in that object's class
- make sure to be able to handle all degenerate nt symbols
- make a generator code that generates many objects for you

*********** 2022-02-24 ************

## install seaborn (for color palettes)

	$ pip3 install seaborn

