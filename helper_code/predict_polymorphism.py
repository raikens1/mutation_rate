#!/usr/bin/env python
from sys import argv
from Bio import SeqIO
from scipy.stats import poisson

usage = """
RA, 6/12/2017
predict_polymorphism
V1
=========================================
Defines a basic predictor object:
Given a fasta file, a dataframe of model parameters, and SNP count
Calculate the poisson binomial parameters and report p count
=========================================
USEAGE: predict_polymorphism.py INCLUDE.bed PARAMS.txt PRIV COUNT(optional)
"""

def main():
    print usage
    if len(argv) == 1:
        print "No arguements specified.  Exiting!"
        exit(0)

    pred = Predictor(argv[2], bool(int(argv[3])))

    pred.getParams(argv[1])

    pred.writeParams(argv[1])

    if len(argv) == 5:
        count = int(argv[4])
        p = pred.ppoibin(count)
    
        print """Using Le Cam's poisson approximation, the probability of 
observing %d or fewer polymorphisms on this sequence is %f.
For a more exact calculation, import the poibin parameters to R and
calculate using the DFT-CF method in package poibin.\n""" % (count, p)


"""
	Predictor Object:
	Counts mutations in a private file by sequence context
	Class variables:
        - prob (dict, str:float) mapping from context to mutation probability
        - chrom (str) name of chromosome for bed line we are currently reading
        - ref_genome_chr (str) fasta sequence of self.chrom from hg19
        - compliments (dict) maps complimentary bases to each other
        - priv (bool) if true, consider private mutation model; else cosmo
        - params (list of floats) running list of params for poibin dist
		- flank (int) number of flanking bases of sequence context to consider
"""
class Predictor(object):
    def __init__(self, modelfile, priv):
        self.prob = {}
        self.chrom = ''
        self.ref_genome_chr = ''
        self.compliments = {"A":"T", "T":"A", "G":"C", "C":"G"}
        self.priv = priv
        self.params = []
        self.mean = 0

        self.readProb(modelfile) # build self.prob
        self.flank = (len(self
            .prob.keys()[0])-1)/2
        print self.prob

    # initialize self.prob and self.flank from prob file
    def readProb(self, modelfile):
        print "Reading in model parameters from %s.\n" % modelfile
        with open(modelfile, 'r') as m:
            m.next() #skip header
            for line in m:
                self.parserow(line.split('\t'))


    # given a row of a modelfile, add to self.prob
    def parserow(self, row):
        if self.priv:
        	self.prob[row[0]] = sum([float(n) for n in row[4:7]])
        else:
            self.prob[row[0]] = sum([float(n) for n in row[1:4]])

    # calculate parameters for sequence
    def getParams(self, bedfile):
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

    # Use Le Cam's theorem to approximate cdf of poibin
    def ppoibin(self, count):
        return poisson.cdf(count, self.mean)

if __name__ == '__main__':
    main()
