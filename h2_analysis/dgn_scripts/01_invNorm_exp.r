####by Heather E. Wheeler 20150302####
date <- Sys.Date()

###############################################
### Directories & Variables
my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
ct.dir <- "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/cis.v.trans.prediction/"
tis <- "DGN-WB"

### Functions
"%&%" = function(a,b) paste(a,b,sep="")
med <- function(x) median(x,na.rm=T)
source(my.dir %&% 'GenABEL/R/ztransform.R')
source(my.dir %&% 'GenABEL/R/rntransform.R')

###############################################
### Directories & Variables
my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
ct.dir <- "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/cis.v.trans.prediction/"
tis <- "DGN-WB"

###Scan data
rpkmid <- ct.dir %&% tis %&% ".exp.ID.list"
expid <- scan(rpkmid,"character")
rpkmgene <- ct.dir %&% tis %&% ".exp.GENE.list"
geneid <- scan(rpkmgene,"character")
rpkmfile <- ct.dir %&% tis %&% ".exp.IDxGENE"
expdata <- scan(rpkmfile)
expdata <- matrix(expdata, ncol = length(geneid), byrow=TRUE)
rownames(expdata) <- expid
colnames(expdata) <- geneid

texpdata <- t(expdata)

##inverse normalize the exp data
rn.expdata<-apply(texpdata,1,"rntransform") 
rn.expdata<-round(rn.expdata,6) #round to 6 places

write.table(rn.expdata,file=my.dir %&% "expArch_DGN-WB_imputedGTs/dgn-exp/" %&% tis %&% ".rntransform.exp.IDxGENE",quote=F,col.names=F,row.names=F)
write(expid,file=my.dir %&% "expArch_DGN-WB_imputedGTs/dgn-exp/" %&% tis %&% "exp.ID.list",ncol=1)
write(geneid,file=my.dir %&% "expArch_DGN-WB_imputedGTs/dgn-exp/" %&% tis %&% "exp.GENE.list",ncol=1)
