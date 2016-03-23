#!/bin/bash
## cat cross-tissue output from 20_GTEx_cross-tissue_CV_elasticNet.r

head -n 1 subsets/cross-tissue_exp_10-foldCV_elasticNet_alpha0.05_hapmapSnpsCEU_chr1_2015-05-15.txt > header
cat subsets/cross-tissue_exp_10-foldCV_elasticNet_alpha0.05_hapmapSnpsCEU_chr* >o
grep -v gene o > p
cat header p > cross-tissue_exp_10-foldCV_elasticNet_alpha0.05_hapmapSnpsCEU_chr1-22_2015-05-15.txt
rm o p header

head -n 1 subsets/cross-tissue_exp_10-foldCV_elasticNet_alpha0.5_hapmapSnpsCEU_chr1_2015-05-14.txt > header
cat subsets/cross-tissue_exp_10-foldCV_elasticNet_alpha0.5_hapmapSnpsCEU_chr* >o
grep -v gene o > p
cat header p > cross-tissue_exp_10-foldCV_elasticNet_alpha0.5_hapmapSnpsCEU_chr1-22_2015-05-14.txt
rm o p header

head -n 1 subsets/cross-tissue_exp_10-foldCV_elasticNet_alpha0.95_hapmapSnpsCEU_chr1_2015-05-15.txt > header
cat subsets/cross-tissue_exp_10-foldCV_elasticNet_alpha0.95_hapmapSnpsCEU_chr* >o
grep -v gene o > p
cat header p > cross-tissue_exp_10-foldCV_elasticNet_alpha0.95_hapmapSnpsCEU_chr1-22_2015-05-15.txt
rm o p header

head -n 1 subsets/cross-tissue_exp_10-foldCV_elasticNet_alpha1_hapmapSnpsCEU_chr1_2015-05-14.txt > header
cat subsets/cross-tissue_exp_10-foldCV_elasticNet_alpha1_hapmapSnpsCEU_chr* >o
grep -v gene o > p
cat header p > cross-tissue_exp_10-foldCV_elasticNet_alpha1_hapmapSnpsCEU_chr1-22_2015-05-14.txt
rm o p header

###pull out R2 for ggplot2

cut -f 1,7 cross-tissue_exp_10-foldCV_elasticNet_alpha0.05_hapmapSnpsCEU_chr1-22_2015-05-15.txt > a
cut -f 7 cross-tissue_exp_10-foldCV_elasticNet_alpha0.5_hapmapSnpsCEU_chr1-22_2015-05-14.txt > b
cut -f 7 cross-tissue_exp_10-foldCV_elasticNet_alpha0.95_hapmapSnpsCEU_chr1-22_2015-05-15.txt > c
cut -f 7 cross-tissue_exp_10-foldCV_elasticNet_alpha1_hapmapSnpsCEU_chr1-22_2015-05-14.txt > d

awk '{print $1,"cross-tissue",$2}' a > e
paste -d' ' e b c d > cross-tissue_exp_10-foldCV_elasticNet_R2_for_ggplot2.txt

rm a b c d e
