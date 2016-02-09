####by Heather E. Wheeler 20141114####
args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()

### hapmap2 non-ambiguous SNPs were extracted from the GTEx genotype data on genegate with:
### /nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/1_vcf2dosage.mach_gtex_hapmapSNPs.pl
### and scp'd to tarbell: /group/im-lab/hwheeler/cross-tissue/gtex-genotypes/

my.dir <- "/group/im-lab/hwheeler/cross-tissue/"
rna.dir <- my.dir %&% "gtex-rnaseq/"
annot.dir <- my.dir %&% "gtex-annot/"
gt.dir <- my.dir %&% "gtex-genotypes/"
grm.dir <- my.dir %&% "gtex-grms/"

gencodefile <- annot.dir %&% "gencode.v18.genes.patched_contigs.summary.protein.chr" %&% args[1]
gencodeset <- args[1]

bimfile <- gt.dir %&% "GTEx_Analysis_2014-06-13.hapmapSnpsCEU.bim" ###get SNP position information###
bim <- read.table(bimfile)
rownames(bim) <- bim$V2

###make globalGRM, only need to do once

machpre <- "GTEx_Analysis_2014-06-13.hapmapSnpsCEU." 

#runGCTAglo <- "gcta64 --dosage-mach-gz " %&% gt.dir %&% machpre %&% "mldose.gz " %&% gt.dir %&% machpre %&% "mlinfo.gz --make-grm-bin --out " %&% grm.dir %&% "GTEx.global"
#system(runGCTAglo)

##make chrGRM

runGCTAchr <- "gcta64 --dosage-mach-gz " %&% gt.dir %&% machpre %&% "chr" %&% gencodeset %&%  ".mldose.gz " %&% gt.dir %&% machpre %&% "chr" %&% gencodeset %&% ".mlinfo.gz --make-grm-bin  --out " %&% grm.dir %&% "GTEx.chr" %&% gencodeset
system(runGCTAchr)

###make localGRMs, for subset of genes in gencodeset

gencode <- read.table(gencodefile)
rownames(gencode) <- gencode[,5]

finished.grms <- scan("done.grms","character") ###already calculated a bunch of grms using old script with whole genome mach files, don't run them again

#ensidlist <- gencode[,5]
ensidlist <- setdiff(rownames(gencode),finished.grms)

for(i in 1:length(ensidlist)){
    cat(i,"/",length(ensidlist),"\n")
    gene <- ensidlist[i]
    geneinfo <- gencode[gene,]
    chr <- geneinfo[1]
    c <- substr(chr$V1,4,5)
    start <- geneinfo$V3 - 1e6 ### 1Mb lower bound for cis-eQTLS
    end <- geneinfo$V4 + 1e6 ### 1Mb upper bound for cis-eQTLs
    chrsnps <- subset(bim,bim[,1]==c) ### pull snps on same chr
    cissnps <- subset(chrsnps,chrsnps[,4]>=start & chrsnps[,4]<=end) ### pull cis-SNP info
    snplist <- cissnps[,2]    
    write.table(snplist, file= my.dir %&% "tmp.SNPlist." %&% gencodeset,quote=F,col.names=F,row.names=F)
    runGCTAgrm <- "gcta64 --dosage-mach-gz " %&% gt.dir %&% machpre %&% "chr" %&% gencodeset %&%  ".mldose.gz " %&% gt.dir %&% machpre %&% "chr" %&% gencodeset %&% ".mlinfo.gz --make-grm-bin --extract tmp.SNPlist." %&% gencodeset %&% " --out " %&% grm.dir %&% gene
    system(runGCTAgrm)
}
