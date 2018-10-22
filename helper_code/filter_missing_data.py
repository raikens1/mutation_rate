#!/usr/bin/env python

import gzip
from sys import argv

"""
RCA, 8/19/2018
filter_missing_data.py
V1
=============================================================
Filters zipped vcf to exclude sites with > 20% missing data
This script was written to filter the SGDP dataset for
replication of the results in 1KG.
USEAGE: python filter_missing_data.py INPUT.vcf.gz
=============================================================
"""


"""
Generate an output vcf.gz which is identical to the input vcf.gz with any rows
with more than 20% missing data removed
"""
def main():
    min_data = 0.8
    infile = argv[1]
    outfile = infile.split(".")[0] + "_missing_filtered.vcf"

    accept = 0
    singles = 0
    all_ref = 0
    missing = 0
    multiallelic = 0


    print "Copying variants with >=%d percent non-missing data from %s to %s" % (min_data*100, infile, outfile)

    o = open(outfile, "a")

    with gzip.open(infile) as f:
	for line in f:
	    if line.startswith('#'):
                o.write(line)
            else:
                row = line.split("\t")
                GTs = parse_row(row)
                if sufficient_data(GTs, min_data):
                    alts = GTs.count("1")
                    if len(row[4]) != 1:
                        multiallelic += 1
                    elif alts == 1:
                        singles += 1
                    elif alts == 0:
                        all_ref +=1
                    else:
                        accept += 1
                        o.write(line)
                else:
                    missing += 1

    o.close()

    reject = multiallelic + missing + singles + all_ref

    print "After filtering, removed %d out of %d variants.  %d remain." % (reject, accept+reject, accept)
    print "\t Multiallelic: ", multiallelic
    print "\t Too much missing data: ", missing
    print "\t No alternate alleles in dataset: ", all_ref
    print "\t Singletons: ", singles

"""
Given a list of alleles, return true if a sufficient amount of data is present
input: GT_list (list), a list of alleles at a locus
input: min_data (float), the percent of data which must be non-missing
input: missing_char (str), character representing missing data
output: (boolean): true if >= min_data % of alleles are not missing
"""
def sufficient_data(GT_list, min_data, missing_char = "."):
    percent_present = 1 - GT_list.count(missing_char)/float(len(GT_list))
    return (percent_present >= min_data)

"""
Parses a vcf line into a list of alleles
input: row (list), elements in a line of a vcf (split by \t characters)
output: (list) a list of the alleles of all individuals in the vcf
"""
def parse_row(row):
    GT_cols = row[9:] #
    GTs = [col[:3] for col in GT_cols] # extract just GT
    allele_lists = [s.split("/") for s in GTs] # separate GTs into lists of alleles
    return [item for sublist in allele_lists for item in sublist]

if __name__ == '__main__':
	main()
