#!/bin/bash

#determine cut-points for SNP chunks required to calculate LDAK LD weights

for i in {1..22}
do
    cp /group/im-lab/nas40t2/Data/Transcriptome/WB1K/imputed/DGN-imputed-for-PrediXcan/DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr${i}.SNPxID dgn-sp-format/DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr${i}.sp
    cp /group/im-lab/nas40t2/Data/Transcriptome/WB1K/imputed/DGN-imputed-for-PrediXcan/DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr${i}.bim dgn-sp-format/DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr${i}.bim
    cp /group/im-lab/nas40t2/Data/Transcriptome/WB1K/imputed/DGN-imputed-for-PrediXcan/DGN.hapmap2.chr1-22.QC.fam dgn-sp-format/DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr${i}.fam

    ./ldak.4.9 --sp dgn-sp-format/DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr${i} --cut-weights chr${i}
done
