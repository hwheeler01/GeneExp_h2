#!/bin/bash

grep -vh --no-filename gene GTEx_subsets/cross-tissue*txt > o
cat header o > cross-tissue_exp_BSLMM-s100K_iterations_all_chr1-22_2015-08-06.txt
    
grep -vh --no-filename gene GTEx_subsets/Adipose\ -\ Subcutaneous*txt > o
cat header o > Adipose-Subcutaneous_TS_exp_BSLMM-s100K_iterations_all_chr1-22_2015-08-06.txt

grep -vh --no-filename gene GTEx_subsets/Artery\ -\ Tibial*txt > o
cat header o > Artery-Tibial_TS_exp_BSLMM-s100K_iterations_all_chr1-22_2015-08-06.txt

grep -vh --no-filename gene GTEx_subsets/Heart\ -\ Left\ Ventricle*txt > o
cat header o > Heart-LeftVentricle_TS_exp_BSLMM-s100K_iterations_all_chr1-22_2015-08-06.txt

grep -vh --no-filename gene GTEx_subsets/Lung_exp*txt > o
cat header o > Lung_TS_exp_BSLMM-s100K_iterations_all_chr1-22_2015-08-06.txt

grep -vh --no-filename gene GTEx_subsets/Muscle\ -\ Skeletal*txt > o
cat header o > Muscle-Skeletal_TS_exp_BSLMM-s100K_iterations_all_chr1-22_2015-08-06.txt

grep -vh --no-filename gene GTEx_subsets/Nerve\ -\ Tibial*txt > o
cat header o > Nerve-Tibial_TS_exp_BSLMM-s100K_iterations_all_chr1-22_2015-08-06.txt

grep -vh --no-filename gene GTEx_subsets/Skin\ -\ Sun\ Exposed*txt > o
cat header o > Skin-SunExposed\(Lowerleg\)_TS_exp_BSLMM-s100K_iterations_all_chr1-22_2015-08-06.txt

grep -vh --no-filename gene GTEx_subsets/Thyroid_exp*txt > o
cat header o > Thyroid_TS_exp_BSLMM-s100K_iterations_all_chr1-22_2015-08-06.txt

grep -vh --no-filename gene GTEx_subsets/Whole\ Blood*txt > o
cat header o > WholeBlood_TS_exp_BSLMM-s100K_iterations_all_chr1-22_2015-08-06.txt
