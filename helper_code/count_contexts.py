#!/usr/bin/env python

from class_counter import Counter
from sys import argv

"""
RA, 7/27/2017
count_contexts.py
V2
===========================================================
Given a SNP list or vcf file and a sequence context window,
Counts the number of occurances of each possible SNP type 
	within that sequence context paradigm.
===========================================================
USEAGE: count_contexts.py INPUTFILE FLANK
"""

def main():

	intro() #print program name and usage

	if len(argv) == 1:
		print "No inputs detected.  Exiting!\n"
		exit(0)
	
	#count infile
	infile = argv[1]
	print "See " + infile[:-3] + "_counts.log for progress on this job." 
	with open(infile[:-3] + "_counts.log", "w+") as log:
		log.write("Tracking progress for this count job.\n")

	counter = Counter(int(argv[2]))

	print counter.counts
	print len(counter.counts)
	counter.count(infile)
	counter.write_counts(infile[:-3] + "_folded_counts.txt")

def intro():
	border =  "=========================================================================\n"
	print border + "\t\t\tCount Contexts!" + "\n"
	print "This tool counts mutations of different types in a given chromosome file.\n"
	print "USEAGE: count_contexts.py INPUTFILE FLANK \n" + border


if __name__ == '__main__':
	main()
