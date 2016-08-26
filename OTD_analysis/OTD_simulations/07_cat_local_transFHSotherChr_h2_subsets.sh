#!/bin/bash
### concatenate output from: 
### 16_est_cross-tissue_local_and_trans_h2.r
### 17_est_tissue-specific_local_and_trans_h2.r
### 18_est_tissue-wide_local_and_trans_h2.r

IFS=$'\n'

## cat tissue-specific
for sim in `cat simlist`;
do
    for tis in `cat nine.spaces.tissue.list`; ##edit tissue.list to add more tissues
    do
	echo $tis
	head -n 1 GTEx.TS.${tis}_${sim}.h2.all.models_FHSfdr0.05.Chr10_globaleQTLOtherChr.2016-0*.txt > header
	cat GTEx.TS.${tis}_${sim}.h2.all.models_FHSfdr0.05.Chr*_globaleQTLOtherChr.2016-0*.txt > o
	grep -v tissue o > p
	sort -t$'\t' -rgk5 p > o
	cat header o > GTEx.TS.${tis}_${sim}.h2.all.models_FHSfdr0.05.Chr1-22_globalOtherChr.2016-06-22.txt
	rm o p header
    done
done

## cat cross-tissue
for sim in `cat simlist`;
do
    head -n 1 cross-tissue_${sim}.h2.all.models_FHSfdr0.05.Chr9_globalOtherChr.2016-0*.txt > header
    cat cross-tissue_${sim}.h2.all.models_FHSfdr0.05.Chr*_globalOtherChr.2016-0*.txt > o
    grep -v gene o > p
    sort -t$'\t' -rgk5 p > o
    cat header o > GTEx.cross-tissue_${sim}.h2.all.models_FHSfdr0.05.Chr1-22_globalOtherChr.2016-06-22.txt 
    rm o p header
done

