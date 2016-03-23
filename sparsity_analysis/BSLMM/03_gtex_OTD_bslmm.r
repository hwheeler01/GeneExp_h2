####by Heather E. Wheeler 20150108####
date <- Sys.Date()
args <- commandArgs(trailingOnly=T)
args <- c('19','C','Skin - Sun Exposed (Lower leg)')
"%&%" = function(a,b) paste(a,b,sep="")
library(dplyr)

###############################################
### Directories & Variables
my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
exp.dir <- my.dir %&% "lmer.fits/"
annot.dir <- my.dir %&% "gtex-annot/"
gt.dir <- my.dir %&% "gtex-genotypes/"
out.dir <- my.dir %&% "BSLMM_exp/GTEx_subsets/"
bslmm.out.dir <- my.dir %&% "BSLMM_exp/"

tis <- args[3]
chrom <- as.numeric(args[1]) 
chrname <- "chr" %&% chrom
whichlist <- args[2]

getquant <- function(x) quantile(x,c(0.5,0.025,0.975)) ##pulls the median and the 95% credible sets

################################################
if(tis == 'cross-tissue'){
  rpkmid <- exp.dir %&% "ranef.mean0.1.SUBJID"
  rpkmgene <- exp.dir %&% "ranef.mean0.1.GENE"
  rpkmfile <- exp.dir %&% "ranef.mean0.1.SUBJIDxGENE"
}else{
  rpkmid <- exp.dir %&% "resid.mean0.1.SAMPID"
  rpkmgene <- exp.dir %&% "resid.mean0.1.GENE"
  rpkmfile <- exp.dir %&% "resid.mean0.1.SAMPIDxGENE"
}
expid <- scan(rpkmid,"character")
geneid <- scan(rpkmgene,"character")
expdata <- scan(rpkmfile) 
expdata <- matrix(expdata, ncol = length(geneid), byrow=TRUE)
rownames(expdata) <- expid
colnames(expdata) <- geneid

t.expdata <- expdata #don't need to transpose

gencodefile <- annot.dir %&% "gencode.v18.genes.patched_contigs.summary.protein"
gencode <- read.table(gencodefile)
rownames(gencode) <- gencode[,5] ##ensid
gencode <- gencode[gencode[,1]==chrname,] ##pull genes on chr of interest
t.expdata <- t.expdata[,intersect(colnames(t.expdata),rownames(gencode))] ###pull gene expression data w/gene info
                
tissues <- read.table (annot.dir %&% "GTEx_Analysis_2014-06-13.SampleTissue.annot",header=T,sep="\t")

if(tis != 'cross-tissue'){
  samplelist <- subset(tissues,SMTSD == tis)
  tissue.exp <- t.expdata[intersect(rownames(t.expdata),samplelist$SAMPID),] ###pull expression data for chosen tissue###

  ###take mean of exp for samples with >1 RNA Seq dataset
  subjidtable <- read.table(annot.dir %&% "GTEx_Data_2014-06-13_Annotations_SubjectSampleMappingDS.txt",
	header=T,sep="\t")
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
  t.expdata <- tissue.exp.substr
}

expsamplelist <- rownames(t.expdata) ###samples with exp data###

bimfile <- gt.dir %&% "GTEx_Analysis_2014-06-13.hapmapSnpsCEU.chr" %&% chrom %&% ".bim" ###get SNP position information###
bim <- read.table(bimfile)
rownames(bim) <- bim$V2
colnames(bim) <- c("chr","snp","cm","bp","ref","alt")

famfile <- gt.dir %&% "GTEx_Analysis_2014-06-13.hapmapSnpsCEU.ID.list"
fam <- scan(famfile,"character")
samplelist <- intersect(fam,expsamplelist)                
                        
exp.w.geno <- t.expdata[samplelist,] ###get expression of samples with genotypes###
explist <- colnames(exp.w.geno)

ngroups = 16 ##split chr into ngroups, edit run*.sh files to match groups

if(whichlist == "A"){
  explist <- explist[1:floor(length(explist)/ngroups)]
}else{
  for(i in 2:(ngroups-1)){
    if(whichlist == LETTERS[i]){
      explist <- explist[(floor(length(explist)/ngroups)*(i-1)+1):(floor(length(explist)/ngroups)*i)]
    }
  }
}
if(whichlist == LETTERS[ngroups]){
  explist <- explist[(floor(length(explist)/ngroups)*(ngroups-1)+1):length(explist)]
}

gtfile <- gt.dir %&% 'GTEx_Analysis_2014-06-13.hapmapSnpsCEU.chr' %&% chrom %&% '.mldose.gz'
gtX <- read.table(gtfile)
a<-gtX[,3:dim(gtX)[2]]
gtX <- as.matrix(a)
colnames(gtX) <- bim$snp
rownames(gtX) <- fam
X <- gtX[samplelist,]
X <- t(X)

resultsarray <- array(0,c(length(explist),19))
dimnames(resultsarray)[[1]] <- explist
resultscol <- c("gene","h50","pve50","rho50","pge50","pi50","n_gamma50","h025","pve025","rho025","pge025",
	"pi025","n_gamma025","h975","pve975","rho975","pge975","pi975","n_gamma975")
dimnames(resultsarray)[[2]] <- resultscol

working100K <- out.dir %&% "working_" %&% tis %&% "_exp_BSLMM-s100K_iterations_chr" %&% chrom %&% whichlist %&% 
	"_" %&% date %&% ".txt"
write(resultscol,file=working100K,ncolumns=19,sep="\t")

for(i in 1:length(explist)){
  cat(i,"/",length(explist),"\n")
  gene <- explist[i]
  geneinfo <- gencode[gene,]
  chr <- geneinfo[1]
  c <- substr(chr$V1,4,5)
  start <- geneinfo$V3 - 1e6 ### 1Mb lower bound for cis-eQTLS
  end <- geneinfo$V4 + 1e6 ### 1Mb upper bound for cis-eQTLs
  chrsnps <- subset(bim,bim[,1]==c) ### pull snps on same chr
  cissnps <- subset(chrsnps,chrsnps[,4]>=start & chrsnps[,4]<=end) ### pull cis-SNP info
  cisgenos <- X[intersect(rownames(X),cissnps[,2]),,drop=FALSE] ### pull cis-SNP genotypes
  if(dim(cisgenos)[1] > 0){
    cisgenos <- mutate(data.frame(cisgenos),snp=rownames(cisgenos))
    cisbim <- inner_join(bim,cisgenos,by='snp')
    annotfile <- cbind(cisbim[,1],cisbim[,4],cisbim[,2])
    genofile <- cbind(cisbim[,1],cisbim[,5:dim(cisbim)[2]])
    phenofile <- data.frame(exp.w.geno[,gene])

    write.table(annotfile, file=out.dir %&% "tmp.annot." %&% chrom %&% ".s." %&% whichlist, quote=F, row.names=F, 
	col.names=F, sep=",")
    write.table(genofile, file=out.dir %&% "tmp.geno." %&% chrom %&% ".s." %&% whichlist, quote=F, row.names=F,
	col.names=F, sep=",")
    write.table(phenofile, file=out.dir %&% "tmp.pheno." %&% chrom %&% ".s." %&% whichlist, quote=F, row.names=F,
	col.names=F, sep=",")

    runBSLMM <- "gemma -g " %&% out.dir %&% "tmp.geno." %&% chrom %&% ".s." %&% whichlist %&% " -p " %&% 
	out.dir %&% "tmp.pheno." %&% chrom %&% ".s." %&% whichlist %&% " -a " %&% out.dir %&% "tmp.annot." %&% 
	chrom %&% ".s." %&% whichlist %&% " -bslmm 1 -seed 12345 -s 100000 -o tmp." %&% chrom %&% ".s." %&% whichlist
    system(runBSLMM)

    hyp <- read.table(bslmm.out.dir %&% "output/tmp." %&% chrom %&% ".s." %&% whichlist %&% ".hyp.txt",header=T)
    hyp50 <- hyp[(dim(hyp)[1]/2+1):dim(hyp)[1],] #take second half of sampling iterations
    quantres <- apply(hyp50,2,getquant)
    res <- c(gene,quantres[1,],quantres[2,],quantres[3,])

  }else{
    res <- c(gene,rep(NA,18))
  }
  names(res) <- c("gene","h50","pve50","rho50","pge50","pi50","n_gamma50","h025","pve025","rho025","pge025",
	"pi025","n_gamma025","h975","pve975","rho975","pge975","pi975","n_gamma975") 
  resultsarray[gene,] <- res
  write(res,file=working100K,ncolumns=19,append=T,sep="\t")
}

write.table(resultsarray,file=out.dir %&% tis %&% "_exp_BSLMM-s100K_iterations_chr" %&% chrom %&% whichlist %&% 
	"_" %&% date %&% ".txt",quote=F,row.names=F,sep="\t")
