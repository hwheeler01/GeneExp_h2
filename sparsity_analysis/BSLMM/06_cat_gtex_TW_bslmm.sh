#!/bin/bash

for TIS in `cat /group/im-lab/nas40t2/hwheeler/cross-tissue/nine.tissue.list`

do

    grep -vh --no-filename gene GTEx_subsets/${TIS}*TW*txt > o
    cat header o > ${TIS}_TW_exp_BSLMM-s100K_iterations_all_chr1-22_2015-10-18.txt


done
