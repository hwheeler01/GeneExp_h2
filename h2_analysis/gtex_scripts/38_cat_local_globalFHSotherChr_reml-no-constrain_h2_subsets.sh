#!/bin/bash
### concatenate output from: 
###35_est_cross-tissue_local_and_globalFHSotherChr_reml-no-constrain_h2.r
###36_est_tissue-specific_local_and_globalFHSotherChr_reml-no-constrain_h2.r
###37_est_tissue-wide_local_and_globalFHSotherChr_reml-no-constrain_h2.r

IFS=$'\n'

## cat tissue-specific
for tis in `cat tissue.list`; ##edit tissue.list to add more tissues
do
    echo $tis
    head -n 1 subsets/GTEx.TS.${tis}.h2.all.models_FHSfdr0.05.Chr1_globaleQTLOtherChr_reml-no-constrain.2015-12-14.txt > header
    cat subsets/GTEx.TS.${tis}.h2.all.models_FHSfdr0.05.Chr*_globaleQTLOtherChr_reml-no-constrain.2015-12-14.txt > o
    grep -v tissue o > p
    sort -t$'\t' -rgk5 p > o
    cat header o > GTEx.TS.${tis}.h2.all.models_FHSfdr0.05.Chr1-22_globaleQTLOtherChr_reml-no-constrain.2015-12-14.txt
    rm o p header
done

## cat tissue-wide
for tis in `cat gtex-rnaseq/ind-tissues-from-nick/GTEx_PrediXmod.tissue.list`;
do
    echo $tis
    head -n 1 subsets/GTEx.TW.${tis}.h2.all.models_FHSfdr0.05.Chr1_globaleQTLOtherChr_reml-no-constrain.2015-12-14.txt > header
    cat subsets/GTEx.TW.${tis}.h2.all.models_FHSfdr0.05.Chr*_globaleQTLOtherChr_reml-no-constrain.2015-12-14.txt > o
    grep -v tissue o > p
    sort -t$'\t' -rgk5 p > o
    cat header o > GTEx.TW.${tis}.h2.all.models_FHSfdr0.05.Chr1-22_globaleQTLOtherChr_reml-no-constrain.2015-12-14.txt
    rm o p header
done

## cat cross-tissue
head -n 1 subsets/cross-tissue.h2.all.models_FHSfdr0.05.Chr1_globaleQTLOtherChr_reml-no-constrain.2015-12-14.txt > header
cat subsets/cross-tissue.h2.all.models_FHSfdr0.05.Chr*_globaleQTLOtherChr_reml-no-constrain.2015-12-14.txt > o
grep -v gene o > p
sort -t$'\t' -rgk5 p > o
cat header o > cross-tissue.h2.all.models_FHSfdr0.05.Chr1-22_globaleQTLOtherChr_reml-no-constrain.2015-12-14.txt 
rm o p header


