#!/bin/sh

for file in *"$1"*.gz; 

do HAND="${file/%.gz/_count}";
bsub -q voight_long -o "$HAND".out -e "$HAND".err python ../../../../../helper_code/count_contexts.py $file 3 3;

done

