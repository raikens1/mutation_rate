#!/usr/bin/env python
from find_context import get_seq_context_variant
from sys import argv
from os import rename
import gzip

"""
RA, 6/12/2017
bin_by_freq
V1
=======================
given a list of SNPs ordered by allele frequency and a bin size m
Sorts the list into bins of size at least m, based on allele frequency
=======================
USEAGE: bin_by_freq.py SNPfile.gz m
"""

# TODO:
# - [ ] add usage printing

def main():
	readSNPs(argv[1], int(argv[2]))


"""helper function for readSNPs
parses a line for allele frequency"""
def get_AF(SNP):
	return float(SNP.split("\t")[-1][:-1].split("=")[1])


"""Reads lines from input file and sorts them into bin files"""
def readSNPs(SNPfile, m):
	b = 1 #bin counter
	handle = SNPfile.split('.')[0]

	print "Saving bin counts and averages to %s_bins.txt"
	log = open(handle+"_bins.txt", "a")
	log.write("\t".join(["#BIN", "j", "AVG_AF", "LOWER_BOUND", "UPPER_BOUNDs\n"]))

	SNPs = gzip.open(SNPfile)
	SNPs.next() #skip header
	SNP = SNPs.readline()
	freq = get_AF(SNP)
	current = freq


	while (True):
		#start new bin
		print "Saving SNPs to bin %d" % b
		o = gzip.open(handle+'_bin'+str(b)+'.gz', "a")
		o.write("#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO\n")

		LB = freq
		sum_AF = 0
		j = 0

		while (j < m): #while bin is not full
			freq = current
			print "Copying over SNPs with frequency %s to this bin" % freq

			#fill bin with all SNPs that have next largest MAF
			while (current == freq):
				sum_AF += current
				o.write(SNP)
				j += 1

				SNP = SNPs.readline()
				if (SNP == ''): #end of file
					print "End of input file.  Closing!"

					#break loop and return
					loginfo = [str(b), str(j), str(sum_AF/j), str(LB), str(freq) +"\n"]
					log.write("\t".join(loginfo))
					log.close()
					o.close()
					SNPs.close()
					return

				else:
					current = get_AF(SNP)

		#bin is full
		UB = freq
		loginfo = [str(b), str(j), str(sum_AF/j), str(LB), str(UB) +"\n"]
		log.write("\t".join(loginfo))
		o.close()
		b+=1

"""
Some code I didn't use but might want later:
==============================================
#write new name to bin file
newname = handle + "_AF_" + str(sum_AF/(j)) + "_j_" + str(j) + ".gz"
rename(handle+'_bin'+str(b)+'.gz', newname)
"""

if __name__ == '__main__':
	main()
