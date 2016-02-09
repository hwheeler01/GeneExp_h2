####by Heather E. Wheeler 20141114####
args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()
library(dplyr)

###make DGN GRMs: local (unique to gene: localGENE.grm), global (all eQTLs from Framingham (FHS): globalChr.grm), global-local (unique to gene:global-GENE_Chr.grm)
###make DGN global GRMs (all SNPs on a chr)
pre.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
my.dir <- pre.dir %&% "expArch_DGN-WB_imputedGTs/"
rna.dir <- my.dir %&% "dgn-exp/"
annot.dir <- pre.dir %&% "gtex-annot/"
#gt.dir <- "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/transfers/PrediXmod/DGN-WB/DGN-imputation/DGN-imputed-for-PrediXcan/"
gt.dir <- "/group/im-lab/nas40t2/Data/Transcriptome/WB1K/imputed/DGN-imputed-for-PrediXcan/"

gencodefile <- annot.dir %&% "gencode.v18.genes.patched_contigs.summary.protein.chr" %&% args[1]
gencodeset <- args[1]

grm.dir <- my.dir %&% "dgn-grms/"

machpre <- "DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU."

bimfile <- gt.dir %&% machpre %&% "chr" %&% gencodeset %&% ".bim" ###get SNP position information###
bim <- read.table(bimfile)
rownames(bim) <- bim$V2

###make globalGRM for chr

runGCTAglo <- "gcta64 --dosage-mach-gz " %&% gt.dir %&% machpre %&% "chr" %&% gencodeset %&% ".mldose.gz " %&% gt.dir %&% machpre %&% "chr" %&% gencodeset %&% ".mlinfo.gz --make-grm-bin --out " %&% grm.dir %&% "DGN.global_Chr" %&% gencodeset
system(runGCTAglo)

