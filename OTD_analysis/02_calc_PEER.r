#before R call run: R CMD INSTALL ~/R_peer_source_1.3.tgz, see run_02_calc_PEER.sh, not needed 6/18/2016
####by Heather E. Wheeler 20141113####
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()
args <- commandArgs(trailingOnly=T)
###############################################
### Directories & Variables
my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
out.dir <- my.dir %&% "paper-reviewer-requests/OTD_simulations/"
rna.dir <- out.dir %&% "sim_exp/"

Nk <- 15 ##number of peer factors to calculate, recommend 25% of sample size, but no more than 100, GTEx included 15 in pilot analyses
simfile <- args[1]

################################################
### Functions & Libraries

library(peer)
library(preprocessCore)
library(data.table)
source(my.dir %&% 'GenABEL/R/ztransform.R')
source(my.dir %&% 'GenABEL/R/rntransform.R')

################################################
exp <- data.frame(fread(rna.dir %&% simfile %&% ".txt"))
expdata <- exp[,3:dim(exp)[2]]
expidlist <- exp$SAMPID
expgenelist <- colnames(expdata)
expmat <- as.matrix(expdata)
expdata <- as.matrix(t(expmat))
colnames(expdata) <- expidlist

###quantile normalize and transform to standard normal expdata matrix, as in GTEx paper###

qn.expdata <- normalize.quantiles(expdata) ##quantile normalize
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
write.table(factors,file=rna.dir %&% simfile %&% ".meanRPKM0.1." %&% Nk %&% ".PEER.factors." %&% date %&% ".txt", quote=F)

finalinfo <- cbind(exp[,1:2],factors)
colnames(finalinfo) <- c("SAMPID","SUBJID","PF1","PF2","PF3","PF4","PF5","PF6","PF7","PF8","PF9","PF10","PF11","PF12","PF13","PF14","PF15")
write.table(finalinfo,file=rna.dir %&% simfile %&% ".meanRPKM0.1.IDinfo.PFs",sep="\t",quote=F,row.names=F)

weights = PEER_getW(model)
precision = PEER_getAlpha(model)
residuals = PEER_getResiduals(model)

pdf(file=rna.dir %&% simfile %&% ".meanRPKM0.1." %&% Nk %&% ".PEER.factors.plotmodel." %&% date %&% ".pdf")
PEER_plotModel(model)
dev.off()


