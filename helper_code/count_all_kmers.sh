#!/bin/sh
# usage: count_all_kmers.sh POP FLANK

for file in *"${1}"*.gz; 

do HAND="${file/%.gz/_count}";
bsub -q voight_normal -o "$HAND".out -e "$HAND".err count_contexts.py $file ${2} ${2};

done

