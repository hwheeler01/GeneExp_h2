#!/bin/bash
### concatenate output from: 
### 16_est_cross-tissue_local_and_trans_h2.r
### 17_est_tissue-specific_local_and_trans_h2.r
### 18_est_tissue-wide_local_and_trans_h2.r

IFS=$'\n'

## cat tissue-specific
for tis in `cat tissue.list`; ##edit tissue.list to add more tissues
do
    echo $tis
    head -n 1 subsets/GTEx.TS.${tis}.h2.all.models_FHSfdr0.05.Chr1_globalOtherChr.2015-10-0*.txt > header
    cat subsets/GTEx.TS.${tis}.h2.all.models_FHSfdr0.05.Chr*_globalOtherChr.2015-10-0*.txt > o
    grep -v tissue o > p
    sort -t$'\t' -rgk5 p > o
    cat header o > GTEx.TS.${tis}.h2.all.models_FHSfdr0.05.Chr1-22_globalOtherChr.2015-10-06.txt
    rm o p header
done

## cat tissue-wide
for tis in `cat gtex-rnaseq/ind-tissues-from-nick/GTEx_PrediXmod.tissue.list`;
do
    echo $tis
    head -n 1 subsets/GTEx.TW.${tis}.h2.all.models_FHSfdr0.05.Chr1_globalOtherChr.2015-10-0*.txt > header
    cat subsets/GTEx.TW.${tis}.h2.all.models_FHSfdr0.05.Chr*_globalOtherChr.2015-10-0*.txt > o
    grep -v tissue o > p
    sort -t$'\t' -rgk5 p > o
    cat header o > GTEx.TW.${tis}.h2.all.models_FHSfdr0.05.Chr1-22_globalOtherChr.2015-10-06.txt
    rm o p header
done

## cat cross-tissue
head -n 1 subsets/cross-tissue.h2.all.models_FHSfdr0.05.Chr10_globalOtherChr.2015-10-05.txt >header
cat subsets/cross-tissue.h2.all.models_FHSfdr0.05.Chr*_globalOtherChr.2015-10-05.txt > o
grep -v gene o > p
sort -t$'\t' -rgk5 p > o
cat header o > GTEx.cross-tissue.h2.all.models_FHSfdr0.05.Chr1-22_globalOtherChr.2015-10-06.txt 
rm o p header


