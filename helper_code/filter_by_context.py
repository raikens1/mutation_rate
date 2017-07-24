#!/usr/bin/env python

from Bio import SeqIO
from sys import argv
import gzip

"""
RA, 6/29/2017
filter_by_context
V2
=======================
Given a vcf file and a mutation type file
this tool prints all the lines with those mutation types
=======================
USEAGE: filter_by_context.py INPUT.GZ MUTATIONS_FILE CHROM
"""

def main():
	#parse args
	infile = argv[1]
	bed = True #TODO: someday I should set this to be user specified
	cf = CFilter(argv[2], argv[3], bed)
	cf.filter(infile)

"""
	CFilter Object:
	Given a SNP list, outputs a file with only the SNPs of a certain type or types
	Class variables:
		- flank (int) number of flanking bases in the sequence context window
		- muts (list) a list of polymorphism types to include
		- mutsfile (str) name of the file where polymorphism types are saved
		- compliments (dict) maps complimentary bases to each other
		- chrom (int) the chromosome being filtered
		- bed (bool) if true, saves the results as a bed file, rather than a SNP list
"""

class CFilter(object):
	def __init__(self, mutsfile, chrom, bed):
		self.compliments = {"A":"T", "T":"A", "G":"C", "C":"G"}
		self.mutsfile = mutsfile
		self.muts = []
		self.readMuts()
		self.flank = (len(self.muts[0])-4)/2
		self.chrom = chrom
		self.bed = bed
		print "Reading reference genome file...",
		self.ref_genome_chr = SeqIO.read('/project/voight_datasets/hg19/chr'+chrom+'.fa', "fasta")
		print "finished!"

	#initialize self.muts
	def readMuts(self):
		with open(self.mutsfile) as f:
			for line in f:
				self.muts.append(line.rstrip())

	#given the position of a SNP on self.chrom, return sequence context in reference genome
	def get_context(self, pos):
		return str(self.ref_genome_chr.seq)[pos-(self.flank+1):pos+self.flank].upper()

	#given a file and a list of mutations, print out a file of just those mutations
	def filter(self, infile):
		outfile = infile[:-3] + "_" + self.mutsfile.split(".")[0]
		print "Printing all polymorphisms of the following types from %s to %s" % (infile, outfile)
		print "\n".join(self.muts)
		
		if self.bed:
			print "Saving results as a .bed file."
			o = open(outfile + ".bed", "a")
			o.write("track name=" + outfile + "\n")
		else:
			print "Saving results as a SNP list file."
			o = gzip.open(outfile, "a")
			o.write("#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO\n")

		with gzip.open(infile) as f:
			for line in f:
				if not line.startswith('#'):
					row = line.split("\t")
					if self.checkrow(row):
						if self.bed:
							line = self.lineToBed(row)
						o.write(line)

		o.close()

	#Converts a line of a SNP list or vcf to simple bed line
	def lineToBed(self,row):
		return "\t".join(["chr" + self.chrom, row[1], str(int(row[1])+1)+"\n"])

	#given a row of a SNP file, return True if that SNP is in self.muts
	def checkrow(self, row):
		context = self.get_context(int(row[1]))
		SNP = context + "->" + row[4]

		if SNP in self.muts:
			return True
		elif self.revc(SNP) in self.muts:
			return True
		else:
			return False

	#given a mutation type, return reverse compliment
	def revc(self, SNP):
		ref = SNP.split("->")[0]
		alt = SNP[-1]

		revc = ''
		for char in ref:
			revc = self.compliments[char] + revc

		return revc + "->" + self.compliments[alt]


if __name__ == '__main__':
	main()
