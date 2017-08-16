#This script merges the SNP lists for CHB, JPT, and CHS into a single SNP list, removing duplicate lines

for number in {1..22}
do

zcat chr"$number"_CHS_EAS.gz chr"$number"_JPT_EAS.gz chr"$number"_CHB_EAS.gz | awk '!seen[$0]++' | bgzip -c > chr"$number"_JPTCHSCHB_EAS.gz

done

zcat chrX_CHS_EAS.gz chrX_JPT_EAS.gz chrX_CHB_EAS.gz | awk '!seen[$0]++' | bgzip -c >> chrX_JPTCHSCHB_EAS.gz
