
# RA, 6/6/2017
# cut_paste_counts.sh
# V1
# =======================
# given a population, 
# reprints outputs of count_contexts.py into a single file which is usable in R
# =======================
# USEAGE: cut_paste_counts.sh POP

paste chr1_"$1"_folded_counts.txt chr2_"$1"_folded_counts.txt  chr3_"$1"_folded_counts.txt  chr4_"$1"_folded_counts.txt  chr5_"$1"_folded_counts.txt  chr6_"$1"_folded_counts.txt  chr7_"$1"_folded_counts.txt  chr8_"$1"_folded_counts.txt  chr9_"$1"_folded_counts.txt  chr10_"$1"_folded_counts.txt  chr11_"$1"_folded_counts.txt  chr12_"$1"_folded_counts.txt  chr13_"$1"_folded_counts.txt  chr14_"$1"_folded_counts.txt  chr15_"$1"_folded_counts.txt  chr16_"$1"_folded_counts.txt  chr17_"$1"_folded_counts.txt  chr18_"$1"_folded_counts.txt  chr19_"$1"_folded_counts.txt  chr20_"$1"_folded_counts.txt  chr21_"$1"_folded_counts.txt  chr22_"$1"_folded_counts.txt chrX_"$1"_folded_counts.txt > temp."$1"
cut -f1-2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65,68,69 temp."$1" > "$1"_chrom_counts.txt
