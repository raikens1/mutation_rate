#!/usr/bin/env python

from sys import argv
from Bio import SeqIO
from itertools import product

usage = """
RA, 8/10/2017
count_gw_contexts.py
V1
=====================================================================
Given a fasta file, a bed file, and a number of flanking nucleotides
Counts the number of contexts within those inclusion regions
	within that sequence context paradigm.
=====================================================================
USEAGE: count_gw_contexts.py BED FLANK
"""
"""
Acknowledgements: This code makes use of some lines helpfully offered
by Onur Yoruk based on his replication timing work!
"""

# TODO:
# - [ ] write main function
# - [ ] debug

def main():
    print usage
    if len(argv) == 1:
        print "No arguements specified.  Exiting!"
        exit(0)
    c = FCount(int(argv[2]))
    len(c.counts)

"""
	FCount Object:
	Counts sequence contexts in a fasta sequence
	Class variables:
        - chrom (str) name of chromosome for bed line we are currently reading
        - ref_genome_chr (str) fasta sequence of self.chrom from hg19
        - compliments (dict) maps complimentary bases to each other
		- flank (int) number of flanking bases of sequence context to consider
"""
class FCount(object):
    def __init__(self, flank):
        self.flank = flank
        self.counts = self.init_counts()
        self.chrom = ''
        self.ref_genome_chr = ''
        self.compliments = {"A":"T", "T":"A", "G":"C", "C":"G"}

    # calculate parameters for sequence
    def count(self, bedfile):
        print "Calculating poisson binomial parameters for regions in %s.\n" % bedfile
        i = 0
        with open(bedfile, 'r') as b:
            for line in b:
                row = line.split('\t')
                if not line.startswith("chr"): # skip comment line
                    continue
                elif row[0][3:] != self.chrom:
                    if i != 0:
                        print "Read %d lines" % i
                    self.switchChrom(row[0][3:])
                    i = 0
                self.parseBedLine(row)
                i += 1
        print "Read %d lines\n" % i
        print "The expected number of polymorphisms on this sequence is %f.\n" % self.mean

    # load new chrom file
    def switchChrom(self, chrom):
        print "Loading new chromosome sequence...",
        self.ref_genome_chr = SeqIO.read('/project/voight_datasets/hg19/chr'+chrom+'.fa', "fasta")
        print "done!"
        self.chrom = chrom

    # given position, get context
    def get_context(self, pos):
        return str(self.ref_genome_chr.seq)[pos-(self.flank+1):pos+self.flank].upper()

    # given a sequence, return reverse compliment
    def reverse_comp(self, sequence):
        sequence_rc = ''
        for char in sequence:
            sequence_rc = self.compliments[char] + sequence_rc

        return sequence_rc

    # read through a line of a bed file and add parameters for each position
    def parseBedLine(self, row):
        start = int(row[1])
        stop = int(row[2])
        n = 0

        for i in range(start+1, stop+1):
            seq = self.get_context(i)

            # skip if N in ref genome
            if "N" in seq:
                n += 1
                continue

            # append substitution prob to parameter list
            elif seq in self.prob:
                p = self.prob[seq]
                self.params.append(p)
                self.mean += p
            else:
                p = self.prob[self.reverse_comp(seq)]
                self.params.append(p)
                self.mean += p
        if n != 0:
            print "Omitted %d contexts with N in reference genome" % n

    # write self.params as a space delimited file (to be read by R as a vector)
    def writeParams(self, bedfile):
        # make output file name
        outfile = bedfile.split('/')[-1].split(".")[0] # remove file extension
        k = 2*self.flank+1
        outfile = outfile + '_' + str(k) + "mer_poibin.params"

        # write to file
        print "Writing these parameters to %s.\n" % outfile
        with open(outfile, 'w+') as o:
            o.write(" ".join([str(f) for f in self.params]))

    def init_counts(self):
        counts = {}
        central_nucs = ['C','A']
        nucs = list('GATC')
        k = 1 + self.flank*2
        for nuc in central_nucs:
            for flanks in itertools.product(nucs,repeat = k-1):
                context = ''.join(flanks[:k/2]) + nuc + ''.join(flanks[k/2:])
                counts[context] = 0
        return  counts

    def only_main_4_bases(self, seq):
        #this function returns True for empty strings
        return all(nuc in list('GATC') for nuc in seq)

    def matches_kmer_length(self, context):
        return (len(context)==self.flank*2+1)

if __name__ == "__main__":
    main()
