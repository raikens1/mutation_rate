#!/usr/bin/env python

from sys import argv
from Bio import SeqIO
from itertools import product

usage = """
RA, 8/10/2017
site_counter.py
V1
=====================================================================
Given a fasta file, a bed file, and a number of flanking nucleotides
Counts the number of contexts within those inclusion regions
	within that sequence context paradigm.
=====================================================================
USEAGE: site_counter.py BED FLANK
"""
"""
Acknowledgements: This code makes use of some lines helpfully offered
by Onur Yoruk based on his replication timing work!
"""


def main():
    print usage
    if len(argv) == 1:
        print "No arguements specified.  Exiting!"
        exit(0)
    c = SCount(int(argv[2]))

    c.count(argv[1])
    outfile = argv[1].split("/")[-1].split(".")[0] + "_site_counts.txt"
    c.write_counts(outfile)

"""
	SCount Object:
	Counts sequence contexts in a fasta sequence
	Class variables:
        - chrom (str) name of chromosome for bed line we are currently reading
        - ref_genome_chr (str) fasta sequence of self.chrom from hg19
        - compliments (dict) maps complimentary bases to each other
		- flank (int) number of flanking bases of sequence context to consider
"""
class SCount(object):
    def __init__(self, flank):
        self.flank = flank
        self.counts = self.init_counts()
        self.chrom = ''
        self.ref_genome_chr = ''
        self.compliments = {"A":"T", "T":"A", "G":"C", "C":"G"}

    # calculate parameters for sequence
    def count(self, bedfile):
    	k = self.flank*2+1
        print "Counting up %d-mer contexts from regions in %s.\n" % (k, bedfile)
        with open(bedfile, 'r') as b:
            for line in b:
                row = line.split('\t')
                if not line.startswith("chr"): # skip comment line
                    continue
                elif row[0][3:] != self.chrom:
                    self.switchChrom(row[0][3:])
                self.parseBedLine(row)

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
            if not self.only_main_4_bases(seq):
                n += 1
                continue

            # append substitution prob to parameter list
            elif seq[self.flank] in 'CA':
            	self.counts[seq] += 1
            else:
                self.counts[self.reverse_comp(seq)] += 1
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
            for flanks in product(nucs,repeat = k-1):
                context = ''.join(flanks[:k/2]) + nuc + ''.join(flanks[k/2:])
                counts[context] = 0
        return  counts

    def only_main_4_bases(self, seq):
        #this function returns True for empty strings
        return all(nuc in list('GATC') for nuc in seq)

    def matches_kmer_length(self, context):
        return (len(context)==self.flank*2+1)

    #sort counts dictionary alphabetically and write it to a file
    def write_counts(self, outfile):
        counts = [value for (key,value) in sorted(self.counts.items())]
        contexts = sorted(self.counts)
        with open(outfile, "w+") as f:
            f.write("Context\tCount\n")
            for i in range(len(contexts)):
                f.write(contexts[i] + "\t" + str(counts[i]) + "\n")



if __name__ == "__main__":
    main()
