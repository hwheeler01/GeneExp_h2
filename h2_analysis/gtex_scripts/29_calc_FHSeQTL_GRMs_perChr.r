####by Heather E. Wheeler 20150323####
args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()
library(dplyr)

### FHS trans-eQTLs (fdr < 0.05),  non-ambiguous SNPs were extracted from the GTEx genotype data with:
### /group/im-lab/nas40t2/hwheeler/cross-tissue/gtex-genotypes/01_vcf2dosage.mach_gtex_FHStrans.pl

### make GTEx GRMs per chromosome using known cis- or trans-eQTLs from FHS at a given threshold:

thresh <- 'fdr0.05'

my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
rna.dir <- my.dir %&% "gtex-rnaseq/"
annot.dir <- my.dir %&% "gtex-annot/"
gt.dir <- my.dir %&% "gtex-genotypes/"
grm.dir <- my.dir %&% "gtex-grms/FHS_eQTLs_" %&% thresh %&% "/"


eqtlfile <- my.dir %&% "expArch_DGN-WB_imputedGTs/Framingham_eqtl-gene_" %&% thresh %&% "_hapmapSnpsCEU.rsIDlist" ##all FHS eQTLs

machpre <- "GTEx_Analysis_2014-06-13.hapmapSnpsCEU.chr" 

###make GTEx chr GRMs using known cis- and trans-eQTLs from FHS

for(chr in 1:22){
   runGCTAgrm <- "gcta64 --dosage-mach-gz " %&% gt.dir %&% machpre %&% chr %&% ".mldose.gz " %&% gt.dir %&% machpre %&% chr %&% ".mlinfo.gz --make-grm-bin --extract " %&% eqtlfile %&% " --out " %&% grm.dir %&% "FHSeQTL-chr" %&% chr %&% " --thread-num 8"
   system(runGCTAgrm)
}
