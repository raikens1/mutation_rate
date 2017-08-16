#!/bin/sh

for file in chr22_MXL_AMR_private.gz; 
do HAND="${file/%.gz/_count}";
bsub -q voight_long -o "$HAND".out -e "$HAND".err python ../../../../helper_code/count_contexts.py $file 1 1;
done
