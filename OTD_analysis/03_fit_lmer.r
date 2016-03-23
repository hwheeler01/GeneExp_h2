####by Heather E. Wheeler 20141114####
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()

###############################################
### Directories & Variables
my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
rna.dir <- my.dir %&% "gtex-rnaseq/"

################################################
### Functions & Libraries
library(lme4)
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

expidinfo <- read.table(rna.dir %&% "GTEx_Analysis_2014-06-13.RNA-seq.meanRPKM0.1.IDinfo.PCs.PFs",sep="\t",header=T)

###quantile normalize and transform to standard normal expdata matrix, as in GTEx paper###

explist <- subset(rowMeans(expdata), rowMeans(expdata)>0.1) ###pull genes with mean expression > 0.1###
explist <- names(explist)
nz.expdata <- expdata[explist,]

rowtable<-function(x) length(table(x))>2 ##function to determine if more than 2 unique exp levels per gene
nonbin<-apply(nz.expdata,1,rowtable) ##apply to matrix
gt2.expdata <- nz.expdata[nonbin,] ##remove genes with <= 2 unique exp levels from matrix

qn.expdata <- normalize.quantiles(gt2.expdata) ##quantile normalize
rn.qn.expdata <- apply(qn.expdata,1,"rntransform") ##rank transform to normality & transposes##

colnames(rn.qn.expdata)<-rownames(gt2.expdata)
rownames(rn.qn.expdata)<-colnames(gt2.expdata)

ranefmat <- matrix(NA,nrow=length(unique(expidinfo$SUBJID)),ncol=dim(rn.qn.expdata)[2]) ##nrow=450 when include PCs1-3 in lmer
residmat <- matrix(NA,nrow=dim(rn.qn.expdata)[1],ncol=dim(rn.qn.expdata)[2])

for(i in 1:dim(rn.qn.expdata)[2]){
    data <- cbind(expidinfo,rn.qn.expdata[,i])
#    fit <- lmer(rn.qn.expdata[,i] ~ (1|SUBJID) + SMTSD + GENDER + PC1 + PC2 + PC3 + PF1 + PF2 + PF3 + PF4 + PF5 + PF6 + PF7 + PF8 + PF9 + PF10 + PF11 + PF12 + PF13 + PF14 + PF15, data)    
    fit <- lmer(rn.qn.expdata[,i] ~ (1|SUBJID) + SMTSD + GENDER + PF1 + PF2 + PF3 + PF4 + PF5 + PF6 + PF7 + PF8 + PF9 + PF10 + PF11 + PF12 + PF13 + PF14 + PF15, data) ##doesn't differ much from fit with PC1-PC3 based on spot check
    fitranef <- ranef(fit)
    fitresid <- resid(fit)
    ranefmat[,i] <- fitranef$SUBJID[,1]
    residmat[,i] <- fitresid
}

rownames(ranefmat) <- rownames(fitranef$SUBJID)
colnames(ranefmat) <- colnames(rn.qn.expdata)
rownames(residmat) <- names(fitresid)
colnames(residmat) <- colnames(rn.qn.expdata)

write.table(ranefmat, file="lmer.fits/ranef.Cross-tissue.exp.pheno_GTEx_Data_2014-06-13.protein-coding.mean0.1_lmer.ranefSUBJID_fixefSMTSD.gender.PFs_" %&% date %&% ".txt",quote=F)
write.table(residmat, file="lmer.fits/resid.Tissue-specific.exp.pheno_GTEx_Data_2014-06-13.protein-coding.mean0.1_lmer.ranefSUBJID_fixefSMTSD.gender.PFs_" %&% date %&% ".txt",quote=F)

write.table(ranefmat, file="lmer.fits/ranef.mean0.1.SUBJIDxGENE",quote=F,row.names=F,col.names=F)
write(colnames(ranefmat),file="lmer.fits/ranef.mean0.1.GENE",ncolumns=1)
write(rownames(ranefmat),file="lmer.fits/ranef.mean0.1.SUBJID",ncolumns=1)

write.table(residmat, file="lmer.fits/resid.mean0.1.SAMPIDxGENE",quote=F,row.names=F,col.names=F)
write(colnames(residmat),file="lmer.fits/resid.mean0.1.GENE",ncolumns=1)
write(rownames(residmat),file="lmer.fits/resid.mean0.1.SAMPID",ncolumns=1)
