#USAGE: ./postprocess_counted_bins.sh POP NO._BINS
#works for outdated output file format

awk '{ a[FNR] = (a[FNR] ? a[FNR] FS : "") $2 } END { for(i=1;i<=FNR;i++) print a[i] }' $(ls -1 *counts.txt)>"$1"_AF_counts_m_"$2".txt
