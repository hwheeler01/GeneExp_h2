####by Heather E. Wheeler 20141114####
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()
args <- commandArgs(trailingOnly=T)

###############################################
### Directories & Variables
my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
rna.dir <- my.dir %&% "gtex-rnaseq/"
sim.dir <- my.dir %&% "paper-reviewer-requests/OTD_simulations/sim_exp/"
#simfile <- "sim_exp_phenotype_errvar-ct_mult-1_seed-123"
simfile <- args[1]

################################################
### Functions & Libraries
library(lme4)
library(preprocessCore)
library(dplyr)
library(data.table)
source(my.dir %&% 'GenABEL/R/ztransform.R')
source(my.dir %&% 'GenABEL/R/rntransform.R')

################################################
exp <- data.frame(fread(sim.dir %&% simfile %&% ".txt"))
expdata <- exp[,3:dim(exp)[2]]
expidlist <- exp$SAMPID
expgenelist <- colnames(expdata)
expmat <- as.matrix(expdata)
expdata <- as.matrix(t(expmat))
colnames(expdata) <- expidlist

###quantile normalize and transform to standard normal expdata matrix, as in GTEx paper###

qn.expdata <- normalize.quantiles(expdata) ##quantile normalize
rn.qn.expdata <- apply(qn.expdata,1,"rntransform") ##rank transform to normality & transposes##

expidinfo <- read.table(rna.dir %&% "GTEx_Analysis_2014-06-13.RNA-seq.meanRPKM0.1.IDinfo.PCs.PFs",sep="\t",header=T) %>% dplyr::select(SAMPID, SMTSD, GENDER)
pfs <- read.table(sim.dir %&% simfile %&% ".meanRPKM0.1.IDinfo.PFs",sep="\t",header=T)
expidinfo <- left_join(pfs,expidinfo,by='SAMPID')

colnames(rn.qn.expdata)<-rownames(expdata)
rownames(rn.qn.expdata)<-colnames(expdata)

ranefmat <- matrix(NA,nrow=length(unique(expidinfo$SUBJID)),ncol=dim(rn.qn.expdata)[2]) ##nrow=450 
residmat <- matrix(NA,nrow=dim(rn.qn.expdata)[1],ncol=dim(rn.qn.expdata)[2])

for(i in 1:dim(rn.qn.expdata)[2]){
    data <- cbind(expidinfo,rn.qn.expdata[,i])
    fit <- lmer(rn.qn.expdata[,i] ~ (1|SUBJID) + SMTSD + GENDER + PF1 + PF2 + PF3 + PF4 + PF5 + PF6 + PF7 + PF8 + PF9 + PF10 + PF11 + PF12 + PF13 + PF14 + PF15, data)
    fitranef <- ranef(fit)
    fitresid <- resid(fit)
    ranefmat[,i] <- fitranef$SUBJID[,1]
    residmat[,i] <- fitresid
}

rownames(ranefmat) <- rownames(fitranef$SUBJID)
colnames(ranefmat) <- colnames(rn.qn.expdata)
rownames(residmat) <- names(fitresid)
colnames(residmat) <- colnames(rn.qn.expdata)

write.table(ranefmat, file="lmer.fits/ranef.Cross-tissue.exp.pheno_" %&% simfile %&% "_lmer.ranefSUBJID_fixefSMTSD.gender.PFs_" %&% date %&% ".txt",quote=F)
write.table(residmat, file="lmer.fits/resid.Tissue-specific.exp.pheno_" %&% simfile %&% "_lmer.ranefSUBJID_fixefSMTSD.gender.PFs_" %&% date %&% ".txt",quote=F)

write.table(ranefmat, file="lmer.fits/ranef_" %&% simfile %&% "_SUBJIDxGENE",quote=F,row.names=F,col.names=F)
write(colnames(ranefmat),file="lmer.fits/ranef_" %&% simfile %&% "_GENE",ncolumns=1)
write(rownames(ranefmat),file="lmer.fits/ranef_" %&% simfile %&% "_SUBJID",ncolumns=1)

write.table(residmat, file="lmer.fits/resid_" %&% simfile %&% "_SAMPIDxGENE",quote=F,row.names=F,col.names=F)
write(colnames(residmat),file="lmer.fits/resid_" %&% simfile %&% "_GENE",ncolumns=1)
write(rownames(residmat),file="lmer.fits/resid_" %&% simfile %&% "_SAMPID",ncolumns=1)
