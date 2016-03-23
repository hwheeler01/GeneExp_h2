####by Heather E. Wheeler 20150108####
date <- Sys.Date()
args <- commandArgs(trailingOnly=T)
#args <- c('19','C','Skin-SunExposed(Lowerleg)')
"%&%" = function(a,b) paste(a,b,sep="")
library(dplyr)

###############################################
### Directories & Variables
my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
exp.dir <- my.dir %&% "gtex-rnaseq/ind-tissues-from-nick/"
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
rpkmid <- exp.dir %&% "GTEx_PrediXmod." %&% tis %&% ".exp.adj.15PEERfactors.3PCs.gender.ID.list"
rpkmgene <- exp.dir %&% "GTEx_PrediXmod." %&% tis %&% ".exp.adj.15PEERfactors.3PCs.gender.GENE.list"
rpkmfile <- exp.dir %&% "GTEx_PrediXmod." %&% tis %&% ".exp.adj.15PEERfactors.3PCs.gender.IDxGENE.RDS"

expid <- scan(rpkmid,"character")
geneid <- scan(rpkmgene,"character")
expdata <- readRDS(rpkmfile) 
##WRONG, transposes and mislabels matrix:expdata <- matrix(expdata, ncol = length(geneid), byrow=TRUE)
##discovered error 20150903, need to rerun all TW

t.expdata <- expdata #don't need to transpose

gencodefile <- annot.dir %&% "gencode.v18.genes.patched_contigs.summary.protein"
gencode <- read.table(gencodefile)
rownames(gencode) <- gencode[,5] ##ensid
gencode <- gencode[gencode[,1]==chrname,] ##pull genes on chr of interest
t.expdata <- t.expdata[,intersect(colnames(t.expdata),rownames(gencode))] ###pull gene expression data w/gene info
                
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
resultscol <- c("gene","h50","pve50","rho50","pge50","pi50","n_gamma50","h025","pve025","rho025","pge025","pi025",
	"n_gamma025","h975","pve975","rho975","pge975","pi975","n_gamma975")
dimnames(resultsarray)[[2]] <- resultscol

working100K <- out.dir %&% "working_" %&% tis %&% "_TW_exp_BSLMM-s100K_iterations_chr" %&% chrom %&% whichlist %&% 
	"_" %&% date %&% ".txt"
write(resultscol,file=working100K,ncolumns=19,sep="\t")

for(i in 47:length(explist)){
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

    write.table(annotfile, file=out.dir %&% "tmppTW.annot." %&% chrom %&% ".s." %&% whichlist, quote=F, row.names=F, 
	col.names=F, sep=",")
    write.table(genofile, file=out.dir %&% "tmppTW.geno." %&% chrom %&% ".s." %&% whichlist, quote=F, row.names=F, 
	col.names=F, sep=",")
    write.table(phenofile, file=out.dir %&% "tmppTW.pheno." %&% chrom %&% ".s." %&% whichlist, quote=F, row.names=F,
	col.names=F, sep=",")

    runBSLMM <- "gemma -g " %&% out.dir %&% "tmppTW.geno." %&% chrom %&% ".s." %&% whichlist %&% " -p " %&% out.dir %&% 
	"tmppTW.pheno." %&% chrom %&% ".s." %&% whichlist %&% " -a " %&% out.dir %&% "tmppTW.annot." %&% chrom %&% ".s." %&% 
	whichlist %&% " -bslmm 1 -seed 12345 -s 100000 -o tmppTW." %&% chrom %&% ".s." %&% whichlist
    system(runBSLMM)

    hyp <- read.table(bslmm.out.dir %&% "output/tmppTW." %&% chrom %&% ".s." %&% whichlist %&% ".hyp.txt",header=T)
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

write.table(resultsarray,file=out.dir %&% tis %&% "_TW_exp_BSLMM-s100K_iterations_chr" %&% chrom %&% whichlist %&% 
	"_" %&% date %&% ".txt",quote=F,row.names=F,sep="\t")
