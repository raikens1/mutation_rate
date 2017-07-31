from sys import argv
from Bio import SeqIO

usage = """
RA, 6/12/2017
predict_polymorphism
V1
=========================================
Defines a basic predictor object:
Given a fasta file, a dataframe of model parameters, and SNP count
Calculate the poisson binomial parameters and report p count
=========================================
USEAGE: predict_polymorphism.py SEQ.fa PARAMS.txt COUNT PRIV
"""

def main():
    if lenth(argv) == 1:
        print usage
        print "No arguements specified.  Exiting!"
        exit(0)

"""
	Predictor Object:
	Counts mutations in a private file by sequence context
	Class variables:
		- flank (int) number of flanking bases of sequence context to consider
		- seq (string) fasta sequence to calculate probabilities over
        - prob (dict, str:float) mapping from context to mutation probability
		- compliments (dict) maps complimentary bases to each other
"""
class Predictor(object):
    def __init__(self, bed, modelfile, priv):
        self.bedfile = bed
        self.prob = {}
        self.flank = 0
        self.chrom = ''
        self.ref_genome_chr = ''
        self.compliments = {"A":"T", "T":"A", "G":"C", "C":"G"}
        self.priv = priv

        self.readProb(params)

    #initialize self.prob and self.flank from prob file
    def readProb(self, modelfile):
        with open(modelfile, 'r') as m:
            m.next() #skip header
            row = m.readline().split('\t')
            self.flank = len(row[0])
            self.parserow(row)
            for line in m:
                self.parserow(line.split('\t')) 

    #given a row of a modelfile, add to self.prob
    def parserow(self, row):
        if self.priv:
            self.param[row[0]] = sum(float(row[4:7]))
        else:
            self.param[row[0]] = sum(float(row[1:3]))

    #calculate parameters for sequence
    def getParams(self):
        pass

    #given a sequence, return reverse compliment
    def reverse_comp(self, sequence):
        sequence_rc = ''
        for char in sequence:
            sequence_rc = self.compliments[char] + sequence_rc

        return sequence_rc

if __name__ == '__main__':
    main()
