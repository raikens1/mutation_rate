#USAGE: ./postprocess_counted_bins.sh FILE_HANDLE NO._BINS
#works on new output file format
#TODO: check that column number is correct

awk '{print $1}' "$1"{1.."$2"}_folded_counts.txt > "$1"_AF_counts.txt
