import gzip
from itertools import product
from Bio import SeqIO

"""
RA, 7/27/2017
class_counter.py
V1
=============================================================
Defines a basic counter object:
Given a vcf and number of flanking nucleotides to include
Writes a file of the counts of polymorphisms for each context
=============================================================
"""

# TODO:
# - [ ] add functionality that lets you choose to go faster with chromosome knowledge
# - [ ] add read-in from reference file instead of initializing from permutations
# - [ ] get rid of One_mer column
# - [ ] handle when more information than GT is provided in GT column

"""
	Counter Object:
	Counts mutations in a private file by sequence context
	Class variables:
		- flank (int) number of flanking bases of sequence context to consider
		- counts (dict) maps contexts to counts (folded)
		- compliments (dict) maps complimentary bases to each other
"""

class Counter(object):
	def __init__(self, flank):
		self.flank = flank
		self.compliments = {"A":"T", "T":"A", "G":"C", "C":"G"}
		self.counts = self.init_counts()
		self.chrom = ''
		self.ref_genome_chr = ''

	#helper function to initialize counts dictionary
	def init_counts(self):
		counts = {}
		length = 2*self.flank + 2
		for combination in product('ACGT', repeat = length):
			bases = ''.join(combination)
			context = bases[:length-1] + '->' + bases[length-1]
			if bases[self.flank] == bases[length-1]: #remove C->C, A->A, mutations, etc.
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
					if context[self.flank] == context[-1]:
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

    # load new chrom file
	def switchChrom(self, chrom):
		print "Loading new chromosome sequence...",
		self.ref_genome_chr = SeqIO.read('/project/voight_datasets/hg19/chr'+chrom+'.fa', "fasta")
		print "done!"
		self.chrom = chrom

	# given position, get context
	def get_context(self, pos):
		return str(self.ref_genome_chr.seq)[pos-(self.flank+1):pos+self.flank].upper()

	#given a line from a file, return sequence context(e.g. "TCC->T")
	def parse_context(self, line):
		row = line.split("\t")

		file_chrom = row[0].split("r")[-1] # if "chr" is prepended to chromosome field, remove it
		if file_chrom != self.chrom:
			self.switchChrom(file_chrom)

		sequence = self.get_context(int(row[1]))
		context = sequence + "->" + row[4]
		return context

	#sort counts dictionary alphabetically and write it to a file
	def write_counts(self, outfile):
		counts = [value for (key,value) in sorted(self.counts.items())]
		contexts = sorted(self.counts)
		with open(outfile, "w+") as f:
			f.write("Context\tCount\tOne_mer\n")
			for i in range(len(contexts)):
				f.write(contexts[i] + "\t" + str(counts[i]) + "\t" + self.one_mer(contexts[i]) + "\n")

	#set all counts to zero
	def reset(self):
		for context in self.counts:
			self.counts[context] = 0

	#given a mutation type, return reverse compliment
	def reverse_comp(self, context):
		window = 2*self.flank + 1
		sequence = context[0:window]
		alt = context[-1]

		sequence_rc = ''
		for char in sequence:
			sequence_rc = self.compliments[char] + sequence_rc

		return sequence_rc + '->' + self.compliments[alt]

	#given a context, find onemer
	def one_mer(self, context):
		ref = context[self.flank]
		if ref in ['T', 'G']:
			return self.compliments[ref] + "->" + self.compliments[context[-1]]
		else:
			return context[self.flank] + '->' + context[-1]
