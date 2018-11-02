
# RA, 6/6/2017
# cut_paste_counts.sh
# V1
# =======================
# given a population, 
# reprints outputs of count_contexts.py into a single file which is usable in R
# =======================
# USEAGE: cut_paste_counts.sh FILE_HANDLE POP

paste chr_1_"$1"_folded_counts.txt chr_2_"$1"_folded_counts.txt  chr_3_"$1"_folded_counts.txt  chr_4_"$1"_folded_counts.txt  chr_5_"$1"_folded_counts.txt  chr_6_"$1"_folded_counts.txt  chr_7_"$1"_folded_counts.txt  chr_8_"$1"_folded_counts.txt  chr_9_"$1"_folded_counts.txt  chr_10_"$1"_folded_counts.txt  chr_11_"$1"_folded_counts.txt  chr_12_"$1"_folded_counts.txt  chr_13_"$1"_folded_counts.txt  chr_14_"$1"_folded_counts.txt  chr_15_"$1"_folded_counts.txt  chr_16_"$1"_folded_counts.txt  chr_17_"$1"_folded_counts.txt  chr_18_"$1"_folded_counts.txt  chr_19_"$1"_folded_counts.txt  chr_20_"$1"_folded_counts.txt  chr_21_"$1"_folded_counts.txt  chr_22_"$1"_folded_counts.txt chr_X_"$1"_folded_counts.txt > temp."$1"
cut -f1-2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65,68,69 temp."$1" > SGDP_"${2}"_chrom_counts.txt
