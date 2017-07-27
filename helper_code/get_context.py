"""
RA, 7/27/2017
get_context.py
V1
============================================================
Defines a function for getting sequence context from a fasta
Arguements:
	chr (string) the chromososme the variant is on
	pos (int) the position of the variatn
	before (int) number of bases of sequence context to left
	after (int) number of bases of sequence context to right
============================================================
"""

def get_context(chr, pos, before, after):
	pos = int(pos) - before - 1 #find start of sequence needed

	# adjust pos to account for newline characters in fasta
	pos += pos/50 - int(pos%50==0)

	ref_file = "/project/voight_datasets/hg19/chr%s.fa" % chr
	hg = open(ref_file, "r")
	hg.readline() # skip header

	# go to position
	hg.seek(pos, 1)

	#read until k characters are added
	k = 1 + before + after
	seq = ''
	while len(seq) !=k:
		seq += hg.read(1).strip()
	hg.close()

	return seq.upper()
