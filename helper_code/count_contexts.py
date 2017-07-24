#!/usr/bin/env python

from find_context import get_seq_context_variant
from sys import argv
import gzip
from itertools import product

"""
RA, 7/24/2017
count_contexts.py
V1
=======================
Given a SNP list or vcf file and a sequence context window,
Counts the number of occurances of each possible SNP type 
	within that sequence context paradign
=======================
USEAGE: count_contexts.py INPUTFILE BEFORE AFTER
"""

# TODO: 
# [] add functionality that lets you choose between biopython and Varun's code
# [] merge before and after into flank
# [] add usage printing

def main():
	border =  "=========================================\n"
	print border + "\t\tCount Contexts!" + "\n"
	print "This tool counts mutations of different types in a given chromosome file.\n" + border
	
	infile = argv[1]
	print "See " + infile[:-3] + "_counts.log for progress on this job." 
	with open(infile[:-3] + "_counts.log", "w+") as log:
		log.write("Tracking progress for this count job.\n")

	counter = Counter(int(argv[2]), int(argv[3]))
	counter.count(infile)
	counter.write_counts(infile[:-3] + "_folded_counts.txt")

"""
	Counter Object:
	Counts mutations in a private file by sequence context
	Class variables:
		- before (int) number of bases before position to include in context window
		- after (int) number of bases after position
		- counts (dict) maps contexts to counts (folded)
		- compliments (dict) maps complimentary bases to each other
"""

class Counter(object):
	def __init__(self, before, after):
		self.before = before
		self.after = after
		self.compliments = {"A":"T", "T":"A", "G":"C", "C":"G"}
		self.counts = self.init_counts()

	#helper function to initialize counts dictionary
	def init_counts(self):
		counts = {}
		length = self.before + self.after + 2
		for combination in product('ACGT', repeat = length):
			bases = ''.join(combination)
			context = bases[:length-1] + '->' + bases[length-1]
			if bases[self.before] == bases[length-1]: #remove C->C, A->A, mutations, etc.
				pass
			elif self.reverse_comp(context) not in counts:
				counts[context] = 0
		print "Count dictionary initialized"
		return counts

	#given a file, count sequence contexts and fill in self.counts
	def count(self, infile):
		print "Counting contexts from input file"
		i = 0
		with gzip.open(infile) as f:
			for line in f:
				if not line.startswith('#'):
					context = self.parse_context(line)
					if context[self.before] == context[-1]:
						with open(infile[:-3] + "_counts.log", "a") as log:
							log.write("Error: reference matches alternate:\n")
							log.write(line)
					elif "N" in context:
						with open(infile[:-3] + "_counts.log", "a") as log:
							log.write("Error: N in reference genome: %s \n" % context)
							log.write(line)							
					elif context in self.counts:
						self.counts[context] += 1
					else:
						self.counts[self.reverse_comp(context)] += 1
					i += 1
					if i % 1000 == 0:
						with open(infile[:-3] + "_counts.log", "a") as log:
							log.write("Counted %s variants\n" % i)


	#given a line from a file, return sequence context(e.g. "TCC->T")
	def parse_context(self, line):
		row = line.split("\t")
		sequence = get_seq_context_variant(row[0], row[1], self.before, self.after)
		context = sequence + "->" + row[4]
		return context

	#sort counts dictionary alphabetically and write it to a file
	def write_counts(self, outfile):
		counts = [value for (key,value) in sorted(self.counts.items())]
		contexts = sorted(self.counts)
		with open(outfile, "w+") as f:# is this syntax right?
			f.write("Context\tCount\tOne_mer\n")
			for i in range(len(contexts)):
				f.write(contexts[i] + "\t" + str(counts[i]) + "\t" + self.one_mer(contexts[i]) + "\n")

	#set all counts to zero
	def reset(self):
		for context in self.counts:
			self.counts[context] = 0

	#given a mutation type, return reverse compliment
	def reverse_comp(self, context):
		window = self.before + self.after + 1
		sequence = context[0:window]
		alt = context[-1]

		sequence_rc = ''
		for char in sequence:
			sequence_rc = self.compliments[char] + sequence_rc

		return sequence_rc + '->' + self.compliments[alt]

	#given a context, find onemer
	def one_mer(self, context):
		ref = context[self.before]
		if ref in ['T', 'G']:
			return self.compliments[ref] + "->" + self.compliments[context[-1]]
		else:
			return context[self.before] + '->' + context[-1]

if __name__ == '__main__':
	main()
