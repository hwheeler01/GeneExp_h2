####by Heather E. Wheeler 20151005####
args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()
library(dplyr)

thresh <- 'fdr0.05' ##threshold for inclusion of trans FHS SNPs

my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
lmer.dir <- my.dir %&% "lmer.fits/"
annot.dir <- my.dir %&% "gtex-annot/"
grm.dir <- my.dir %&% "gtex-grms/"
grm.trans.dir <- grm.dir %&% "FHS_eQTLs_" %&% thresh %&% "/"

gencodefile <- annot.dir %&% "gencode.v18.genes.patched_contigs.summary.protein.chr" %&% args[1] ##genes split into 22 chrs
gencodeset <- args[1]

###############################################
### Scan expression data
idfile <- lmer.dir %&% "resid.mean0.1.SAMPID"
genefile <- lmer.dir %&% "resid.mean0.1.GENE"
expfile <- lmer.dir %&% "resid.mean0.1.SAMPIDxGENE"

sampid <- scan(idfile,"character")
geneid <- scan(genefile, "character")
expdata <- scan(expfile)
expdata <- matrix(expdata, ncol = length(geneid), byrow=TRUE)
rownames(expdata) <- sampid
colnames(expdata) <- geneid

### Get gene subset to analyze
localfile <- grm.dir %&% "localGRM.list"
locallist <- scan(localfile,"character")

gencode <- read.table(gencodefile)
rownames(gencode) <- gencode[,5]

grmlist <- intersect(locallist,rownames(gencode)) ##genes in genecodeset with grm
ensidlist <- intersect(grmlist,colnames(expdata)) ##genes with grm and exp

gensetexp <- expdata[,ensidlist] ##get exp data from intersected genes

### Get tissue sets to analyze
tissues <- read.table (annot.dir %&% "GTEx_Analysis_2014-06-13.SampleTissue.annot",header=T,sep="\t")
tissuesets <- table(tissues$SMTSD)
largesets <- tissuesets[tissuesets >= 100] ###start with tissues with n>=100
tislist <- names(largesets)

for(h in 1:length(tislist)){
  tis <- tislist[h]
  samplelist <- subset(tissues,SMTSD == tis)
  tissue.exp <- gensetexp[intersect(rownames(gensetexp),samplelist$SAMPID),] ###pull expression data for chosen tissue###
    
  ###take mean of exp for samples with >1 RNA Seq dataset
  subjidtable <- read.table(annot.dir %&% "GTEx_Data_2014-06-13_Annotations_SubjectSampleMappingDS.txt",header=T,sep="\t")
  rownames(subjidtable) <- subjidtable$SAMPID
  subjidlist <- subjidtable[rownames(tissue.exp),1]
  tissue.exp.substr <- tissue.exp
  rownames(tissue.exp.substr) <- subjidlist
  uniqsubjid <- unique(subjidlist)
  for(id in uniqsubjid){  ####take mean of exp for samples with >1 RNA Seq dataset
    matchexp <- tissue.exp.substr[rownames(tissue.exp.substr)==id,]
    if(is.array(matchexp)=='TRUE'){
      expmean <- colMeans(matchexp)
		  tissue.exp.substr[as.character(id),] <- expmean
    }
  }

  ### Get subject subset to analyze
  grm.id <- read.table(grm.dir %&% "GTEx.global.grm.id")
  indidlist <- intersect(uniqsubjid,grm.id[,1]) ##subjects with exp and grm
  nsubj <- length(indidlist)

  ### Get expression data from intersected genes and subjects
  localexp <- tissue.exp.substr[indidlist,]
  localensid <- colnames(localexp) #genelist to analyze

  ### Output matrices
  loc.mat <- matrix(0,nrow=length(localensid),ncol=7)
  colnames(loc.mat) <- c("tissue","N","ensid","gene","local.h2","local.se","local.p")
  
  glo.mat <- matrix(0,nrow=length(localensid),ncol=4)
  colnames(glo.mat) <- c("gene","global.h2","global.se","global.p")
  
  locglo.mat <- matrix(0,nrow=length(localensid),ncol=5)
  colnames(locglo.mat) <- c("gene","loc.jt.h2","loc.jt.se","glo.jt.h2","glo.jt.se")
  
  for(i in 1:length(localensid)){
    cat(i,"of",length(localensid),"\n")
    ensid <- as.character(localensid[i])
    gene <- as.character(gencode[ensid,6])
    chr <- as.character(gencode[ensid,1])
    c <- as.numeric(substr(chr,4,5))
    
    #output expression pheno for gcta
    geneexp <- cbind(rownames(localexp),localexp[,ensid])
    write.table(geneexp, file="tmp.TS.pheno." %&% thresh %&% gencodeset, col.names=F, quote=F) #output pheno for gcta
    
    ## Y ~ localGRM
    runLOC <- "gcta64 --grm " %&% grm.dir %&% ensid %&% " --reml-no-constrain --pheno tmp.TS.pheno." %&% thresh %&% gencodeset %&% " --out tmp.TS." %&% thresh %&% gencodeset
    system(runLOC)
    if(file.exists("tmp.TS." %&% thresh %&% gencodeset %&% ".hsq")==TRUE){
        hsq <- scan("tmp.TS." %&% thresh %&% gencodeset %&% ".hsq","character")
        res <- c(tis, nsubj, ensid, gene, hsq[14], hsq[15], hsq[25])
    }else{
	res <- c(tis, nsubj, ensid, gene, NA, NA, NA) #gcta did not coverge or Error: the information matrix is not invertible.
    }
    loc.mat[i,] <- res
    system("rm tmp.TS." %&% thresh %&% gencodeset %&% ".hsq")
    
    ## Y ~ globalGRM (all FHS eQTLs, not including localGRM SNPs) ##changed 3/18/15 to all FHS eQTLs on other chrs
    ###combine GRMs
    chrlist <- c(1:22)
    otherchrs <- setdiff(chrlist,c)
    grmmat <- matrix(NA,nrow=length(otherchrs),ncol=1)
    for(j in 1:length(otherchrs)){
      chrom <- otherchrs[j]
      grmfile <- grm.trans.dir %&% "FHSeQTL-chr" %&% chrom
      grmmat[j,] <- grmfile
    }
    #	glochrfile <- grm.dir %&% "global-" %&% gene %&% "-Chr" %&% gencodeset
    #	grmmat[length(chrlist),] <- glochrfile
    write.table(grmmat,file="tmp.TS.grm.list" %&% thresh %&% gencodeset,quote=F,col.names=F,row.names=F)
    runGRM <- "gcta64 --mgrm-bin tmp.TS.grm.list" %&% thresh %&% gencodeset %&% " --make-grm --out " %&% grm.trans.dir %&% "global-" %&% ensid %&% "-otherChr"
    system(runGRM)
    
    runGLO <- "gcta64 --grm " %&% grm.trans.dir %&% "global-" %&% ensid %&% "-otherChr --reml-no-constrain --pheno tmp.TS.pheno." %&% thresh %&% gencodeset %&% " --out tmp.TS." %&% thresh %&% gencodeset
    system(runGLO)
    if(file.exists("tmp.TS." %&% thresh %&% gencodeset %&% ".hsq")==TRUE){
	hsq <- scan("tmp.TS." %&% thresh %&% gencodeset %&% ".hsq","character")
        res <- c(gene, hsq[14], hsq[15], hsq[25])
    }else{
	res <- c(gene, NA, NA, NA) #gcta did not coverge or Error: the information matrix is not invertible.
    }
    glo.mat[i,] <- res 
    system("rm tmp.TS." %&% thresh %&% gencodeset %&% ".hsq")
    
    ## Y ~ localGRM + globalGRM (joint model)
    runMULT <- "echo " %&% grm.dir %&% ensid %&% " > tmp.TS.multiGRM." %&% thresh %&% gencodeset
    runMULT2 <- "echo " %&% grm.trans.dir %&% "global-" %&% ensid %&% "-otherChr >> tmp.TS.multiGRM." %&% thresh %&% gencodeset
    system(runMULT)
    system(runMULT2)
    runLOCGLO <- "gcta64 --mgrm-bin tmp.TS.multiGRM." %&% thresh %&% gencodeset %&% " --reml-no-constrain --pheno tmp.TS.pheno." %&% thresh %&% gencodeset %&% " --out tmp.TS." %&% thresh %&% gencodeset
    system(runLOCGLO)
    if(file.exists("tmp.TS." %&% thresh %&% gencodeset %&% ".hsq")==TRUE){
        hsq <- scan("tmp.TS." %&% thresh %&% gencodeset %&% ".hsq","character")
        res <- c(gene, hsq[17], hsq[18], hsq[20], hsq[21])
    }else{
	res <- c(gene, NA, NA, NA, NA) #gcta did not coverge or Error: the information matrix is not invertible.
    }
    locglo.mat[i,] <- res
    system("rm tmp.TS." %&% thresh %&% gencodeset %&% ".hsq")
  }

  full.mat <- cbind(loc.mat,glo.mat,locglo.mat)
  write.table(full.mat,file="GTEx.TS." %&% tis %&% ".h2.all.models_FHS" %&% thresh %&% ".Chr" %&% gencodeset %&% "_globaleQTLOtherChr_reml-no-constrain." %&% date %&% ".txt",quote=F,row.names=F,sep="\t")

}

