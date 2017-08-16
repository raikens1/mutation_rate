#USAGE: ./cat_and_sort.sh POP

zcat *"$1"*.gz | bgzip -cd | cut -d ";" -f1 | sort -k 8.7 | tail -n +23 |bgzip -c > "$1"_sortedSNPs.gz
