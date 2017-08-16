for file in ../monoallelic_variants/*.gz;
    do CHR="${file//[!0-9]/}";
    if [ -z "$CHR" ]; then CHR="X"; fi;
    STR="${file/%.vcf.gz/_private}"
    HAND="${STR:24}"
    vcf-isec -c "$file" ../../../cosmo/cosmo_multiallelic/chr"$CHR"_COSMO_private_multi.vcf.gz -f | cut -f1-9 | bgzip -c > "$HAND".gz;
    done
