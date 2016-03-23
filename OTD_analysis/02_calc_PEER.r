#before R call run: R CMD INSTALL ~/R_peer_source_1.3.tgz, see run_02_calc_PEER.sh
####by Heather E. Wheeler 20141113####
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()

###############################################
### Directories & Variables
my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
rna.dir <- my.dir %&% "gtex-rnaseq/"

Nk <- 15 ##number of peer factors to calculate, recommend 25% of sample size, but no more than 100, GTEx included 15 in pilot analyses

################################################
### Functions & Libraries

library(peer)
library(preprocessCore)

source(my.dir %&% 'GenABEL/R/ztransform.R')
source(my.dir %&% 'GenABEL/R/rntransform.R')

################################################
expidlist <- scan(rna.dir %&%  "GTEx_Analysis_2014-06-13.RNA-seq.ID.list","character")
expgenelist <- scan(rna.dir %&% "GTEx_Analysis_2014-06-13.RNA-seq.GENE.list","character")
exp <- scan(rna.dir %&% "GTEx_Analysis_2014-06-13.RNA-seq.GENExID")
expdata <- matrix(exp, ncol=length(expidlist), byrow=T)
colnames(expdata) <- expidlist
rownames(expdata) <- expgenelist

expidinfo <- read.table(rna.dir %&% "GTEx_Analysis_2014-06-13.RNA-seq.IDinfo.list",sep="\t")

###quantile normalize and transform to standard normal expdata matrix, as in GTEx paper###

explist <- subset(rowMeans(expdata), rowMeans(expdata)>0.1) ###pull genes with mean expression > 0.1###
explist <- names(explist)
nz.expdata <- expdata[explist,]

rowtable<-function(x) length(table(x))>2 ##function to determine if more than 2 exp levels per gene
nonbin<-apply(nz.expdata,1,rowtable) ##apply to matrix
gt2.expdata <- nz.expdata[nonbin,] ##remove genes with <= 2 exp levels from matrix

qn.expdata <- normalize.quantiles(gt2.expdata) ##quantile normalize
rn.qn.expdata <- apply(qn.expdata,1,"rntransform") ##rank transform to normality & transposes##

###Now we can create the model object, ### from https://github.com/PMBio/peer/wiki/Tutorial

model = PEER()

###set the observed data,

PEER_setPhenoMean(model,as.matrix(rn.qn.expdata))

dim(PEER_getPhenoMean(model))

###(NULL response means no error here), say we want to infer K=20 hidden confounders,

PEER_setNk(model,Nk)

PEER_getNk(model)

####and perform the inference. ###for Nk=20 and GTEx-NT, it took 323 iterations, for Nk=15 and GTEx-NT, it took 37 iterations

PEER_update(model)

factors = PEER_getX(model)
rownames(factors) <- colnames(expdata)
write.table(factors,file=rna.dir %&% "All_tissues.GTEx_Data_2014-06-13.protein-coding.meanRPKM0.1." %&% Nk %&% ".PEER.factors." %&% date %&% ".txt", quote=F)

finalinfo <- cbind(expidinfo,factors)
colnames(finalinfo) <- c("SAMPID","SUBJID","SMTSD","GENDER","PC1","PC2","PC3","PF1","PF2","PF3","PF4","PF5","PF6","PF7","PF8","PF9","PF10","PF11","PF12","PF13","PF14","PF15")
write.table(finalinfo,file=rna.dir %&% "GTEx_Analysis_2014-06-13.RNA-seq.meanRPKM0.1.IDinfo.PCs.PFs",sep="\t",quote=F,row.names=F)

weights = PEER_getW(model)
precision = PEER_getAlpha(model)
residuals = PEER_getResiduals(model)

pdf(file=rna.dir %&% "All_tissues.GTEx_Data_2014-06-13.protein-coding.meanRPKM0.1." %&% Nk %&% ".PEER.factors.plotmodel." %&% date %&% ".pdf")
PEER_plotModel(model)
dev.off()


