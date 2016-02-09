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
exp.dir <- my.dir %&% "gtex-rnaseq/ind-tissues-from-nick/" ##tissue exp adj by 15 PEER, 3PCs, gender

gencodefile <- annot.dir %&% "gencode.v18.genes.patched_contigs.summary.protein.chr" %&% args[1] ##genes split into 22 chrs
gencodeset <- args[1]

###############################################

### Get gene subset to analyze
localfile <- grm.dir %&% "localGRM.list"
locallist <- scan(localfile,"character")

gencode <- read.table(gencodefile)
rownames(gencode) <- gencode[,5]

grmlist <- intersect(locallist,rownames(gencode)) ##genes in genecodeset with grm

### Get tissue sets to analyze
tissues <- read.table (annot.dir %&% "GTEx_Analysis_2014-06-13.SampleTissue.annot",header=T,sep="\t")
tissuesets <- table(tissues$SMTSD)
largesets <- tissuesets[tissuesets >= 100] ###start with tissues with n>=100
tislist <- names(largesets)

tiswidefile <- exp.dir %&% "GTEx_PrediXmod.tissue.list"
peerlist<- scan(tiswidefile,"character")
squishlist<- gsub(' ','',tislist) ##removes all whitespace to match peerlist & .RDS files

tislist<-intersect(peerlist,squishlist)

for(h in 1:length(tislist)){
  tis <- tislist[h]

  ##read exp data
  obs<-readRDS(exp.dir %&% "GTEx_PrediXmod." %&% tis %&% ".exp.adj.15PEERfactors.3PCs.gender.IDxGENE.RDS")
  genefile <- exp.dir %&% "GTEx_PrediXmod." %&% tis %&% ".exp.adj.15PEERfactors.3PCs.gender.GENE.list"
  idfile <- exp.dir %&% "GTEx_PrediXmod." %&% tis %&% ".exp.adj.15PEERfactors.3PCs.gender.ID.list"
  genelist <- scan(genefile,"character")
  idlist <- scan(idfile,"character")
  colnames(obs) <- genelist
  rownames(obs) <- idlist    
  tissue.exp <- obs[,intersect(colnames(obs),grmlist)] ###pull expression data for chosen tissue & grmlist###
    
  nsubj <- dim(tissue.exp)[1]
  ### Get expression data from intersected genes and subjects
  localexp <- tissue.exp
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
    write.table(geneexp, file="tmp.TW.pheno." %&% thresh %&% gencodeset, col.names=F, quote=F) #output pheno for gcta
    
    ## Y ~ localGRM
    runLOC <- "gcta64 --grm " %&% grm.dir %&% ensid %&% " --reml --pheno tmp.TW.pheno." %&% thresh %&% gencodeset %&% " --out tmp.TW." %&% thresh %&% gencodeset
    system(runLOC)
    hsq <- scan("tmp.TW." %&% thresh %&% gencodeset %&% ".hsq","character")
    if(hsq[23]=="of"){
      res <- c(tis, nsubj, ensid, gene, NA, NA, NA) ##gcta did not converge, script is reading in previous Y~localGRM+globalGRM result
    }else{
      res <- c(tis, nsubj, ensid, gene, hsq[14], hsq[15], hsq[25])
    }
    loc.mat[i,] <- res
    
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
    write.table(grmmat,file="tmp.TW.grm.list" %&% thresh %&% gencodeset,quote=F,col.names=F,row.names=F)
    runGRM <- "gcta64 --mgrm-bin tmp.TW.grm.list" %&% thresh %&% gencodeset %&% " --make-grm --out " %&% grm.trans.dir %&% "global-" %&% ensid %&% "-otherChr"
    system(runGRM)
    
    runGLO <- "gcta64 --grm " %&% grm.trans.dir %&% "global-" %&% ensid %&% "-otherChr --reml --pheno tmp.TW.pheno." %&% thresh %&% gencodeset %&% " --out tmp.TW." %&% thresh %&% gencodeset
    system(runGLO)
    hsq <- scan("tmp.TW." %&% thresh %&% gencodeset %&% ".hsq","character")
    res <- c(gene, hsq[14], hsq[15], hsq[25])
    glo.mat[i,] <- res 
    
    
    ## Y ~ localGRM + globalGRM (joint model)
    runMULT <- "echo " %&% grm.dir %&% ensid %&% " > tmp.TW.multiGRM." %&% thresh %&% gencodeset
    runMULT2 <- "echo " %&% grm.trans.dir %&% "global-" %&% ensid %&% "-otherChr >> tmp.TW.multiGRM." %&% thresh %&% gencodeset
    system(runMULT)
    system(runMULT2)
    runLOCGLO <- "gcta64 --mgrm-bin tmp.TW.multiGRM." %&% thresh %&% gencodeset %&% " --reml --pheno tmp.TW.pheno." %&% thresh %&% gencodeset %&% " --out tmp.TW." %&% thresh %&% gencodeset
    system(runLOCGLO)
    hsq <- scan("tmp.TW." %&% thresh %&% gencodeset %&% ".hsq","character")
    if(hsq[18]=="logL0"){
      res <- c(gene, NA, NA, NA, NA) ##gcta did not converge, script is reading in Y~globalGRM result
    }else{
      res <- c(gene, hsq[17], hsq[18], hsq[20], hsq[21])
    }
    locglo.mat[i,] <- res
  }
  
  full.mat <- cbind(loc.mat,glo.mat,locglo.mat)
  write.table(full.mat,file="GTEx.TW." %&% tis %&% ".h2.all.models_FHS" %&% thresh %&% ".Chr" %&% gencodeset %&% "_globalOtherChr." %&% date %&% ".txt",quote=F,row.names=F,sep="\t")
  
}

