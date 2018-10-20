To process from step 2 to step 3 (i.e. removing rows with >20% missing data), the following command was run:

`for file in chr_*.vcf.gz; do bsub -q voight_normal -o filter_${file}.out -e filter_${file}.err "filter_missing_data.py ${file}"; done`
