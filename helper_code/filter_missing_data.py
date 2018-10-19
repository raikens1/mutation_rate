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

# TODO:
# - [ ] implement basic functionality
# - [ ] better naming for output file
# - [ ] test

def main():
    # open output file in write mode
    # with input file in read mode:
    #   rewrite header to new file
    #   parse GT information from row
    #   if >=20% of GTs present, print line to new file
    min_data = 0.2
    infile = argv[1]
	outfile = "outfile.vcf.gz"# infile[:-3] + "_" + self.mutsfile.split(".")[0]

	print "Copying variants with <=%d\% missing data from %s to %s" % (min_data*100, infile, outfile)

    o = gzip.open(outfile, "a")

	with gzip.open(infile) as f:
		for line in f:
			if line.startswith('#'):
                o.write(line)
            else:
				row = line.split("\t")
				if sufficient_data(row, min_data):
					o.write(line)

		o.close()
    # bonus points: tabix new file.

def sufficient_data(row, min_data):
    return False


if __name__ == '__main__':
	main()
