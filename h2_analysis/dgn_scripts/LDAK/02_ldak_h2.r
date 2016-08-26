####by Heather E. Wheeler 20141114####
args <- commandArgs(trailingOnly=T)
#args <- "22"
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()
library(dplyr)

pre.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
my.dir <- pre.dir %&% "expArch_DGN-WB_imputedGTs/"
rna.dir <- my.dir %&% "dgn-exp/"
annot.dir <- pre.dir %&% "gtex-annot/"
gt.dir <- "/group/im-lab/nas40t2/Data/Transcriptome/WB1K/imputed/DGN-imputed-for-PrediXcan/"

gencodefile <- annot.dir %&% "gencode.v18.genes.patched_contigs.summary.protein.chr" %&% args[1]
gencodeset <- args[1]

grm.dir <- my.dir %&% "dgn-grms/"

machpre <- "DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU."

bimfile <- gt.dir %&% machpre %&% "chr" %&% gencodeset %&% ".bim" ###get SNP position information###
bim <- read.table(bimfile)
rownames(bim) <- bim$V2

###get dosage data
gt <- read.table(gt.dir %&% "DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr" %&% args[1] %&% ".SNPxID")

### Scan expression data
idfile <- rna.dir %&% "DGN-WBexp.ID.list"
genef <- rna.dir %&% "DGN-WBexp.GENE.list"
expfile <- rna.dir %&% "DGN-WB.rntransform.exp.IDxGENE"

subjid <- scan(idfile,"character")
geneid <- scan(genef, "character")
expdata <- scan(expfile)
expdata <- matrix(expdata, ncol = length(geneid), byrow=TRUE)
rownames(expdata) <- subjid
colnames(expdata) <- geneid

gencode <- read.table(gencodefile)
colnames(gencode) <- c('chr','str','start','end','ensid','gene','func','known')

### get chr gene exp data
genebool <- colnames(expdata) %in% gencode$gene
chrexp <- expdata[,genebool==TRUE]
expgencode <- gencode[gencode$gene %in% colnames(expdata),]
rownames(expgencode) = expgencode$ensid
ensidlist <- expgencode[,5]

#for(i in 1:10){
for(i in 1:length(ensidlist)){
    cat(i,"/",length(ensidlist),"\n")
    ensid <- ensidlist[i]
    geneinfo <- expgencode[as.character(ensid),]
    gene <- geneinfo[1,6]
    chr <- geneinfo[1,1]
    c <- substr(chr,4,5)
    start <- geneinfo$start - 1e6 ### 1Mb lower bound for cis-eQTLS
    end <- geneinfo$end + 1e6 ### 1Mb upper bound for cis-eQTLs
    genefile <- data.frame(as.character(ensid),c,start,end)
    write.table(genefile,file="genefile" %&% c,quote=F,row.names=F,col.names=F)

    exppheno <- data.frame(chrexp[,as.character(gene)])
    colnames(exppheno) <- c('exp')
    phenofile <- mutate(exppheno, FID=rownames(exppheno), IID=1) %>% dplyr::select(FID, IID, exp)
    write.table(phenofile,file="phenofile" %&% c,quote=F,row.names=F)

    ldak1 <- "./ldak.4.9 --cut-genes gene" %&% c %&% " --sp dgn-sp-format/DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr" %&% c %&% " --genefile genefile" %&% c %&% " --weights chr" %&% c %&% "/weightsALL"
    ldak2 <- "./ldak.4.9 --calc-genes-reml gene" %&% c %&% " --pheno phenofile" %&% c %&% " --sp dgn-sp-format/DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr" %&% c %&% " --partition 1 --weights chr" %&% c %&% "/weightsALL"
    system(ldak1)
    system(ldak2)

    res <- read.table("gene" %&% c %&% "/regress1",header=T)
    res <- mutate(res,gene=gene)

    if(exists("allres") == FALSE){
      allres = res
    }else{
      allres <- rbind(allres, res)
    }
}

write.table(allres, file="DGN_ldak_reml_chr" %&% c %&% "_" %&% date %&% ".txt", quote=F, row.names=F)


