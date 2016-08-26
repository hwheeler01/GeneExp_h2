####by Heather E. Wheeler 20160617####
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()
library(dplyr)

my.dir <- "/Volumes/im-lab/nas40t2/hwheeler/cross-tissue/"
lmer.dir <- my.dir %&% "lmer.fits/"
annot.dir <- my.dir %&% "gtex-annot/"
otd.dir <- my.dir %&% "paper-reviewer-requests/OTD_simulations/"

###############################################
### get sample info
infofile <- my.dir %&% "gtex-rnaseq/GTEx_Analysis_2014-06-13.RNA-seq.IDinfo.list"
info <- read.table(infofile,sep="\t")
info <- dplyr::rename(info,SAMPID=V1,SUBJID=V2,TISSUE=V3)

### Scan CT expression data
idfile <- lmer.dir %&% "ranef.mean0.1.SUBJID"
genefile <- lmer.dir %&% "ranef.mean0.1.GENE"
expfile <- lmer.dir %&% "ranef.mean0.1.SUBJIDxGENE"

subjid <- scan(idfile,"character")
geneid <- scan(genefile, "character")
expdata <- scan(expfile)
expdata <- matrix(expdata, ncol = length(geneid), byrow=TRUE)
rownames(expdata) <- subjid
colnames(expdata) <- geneid

### Get individual subset to analyze
grm.id <- read.table(my.dir %&% "gtex-grms/GTEx.global.grm.id")
indidlist <- intersect(subjid,grm.id[,1]) ##subjects with exp and grm
nsubj <- length(indidlist)
tis <- "cross-tissue"

### Get expression data from intersected subjects
ctexp <- expdata[indidlist,]
##output for future use
write.table(ctexp,otd.dir %&% "obs_otd_exp/CT_exp.txt",quote=F)

###############################################
### Scan TS expression data
tsidfile <- lmer.dir %&% "resid.mean0.1.SAMPID"
tsgenefile <- lmer.dir %&% "resid.mean0.1.GENE"
tsexpfile <- lmer.dir %&% "resid.mean0.1.SAMPIDxGENE"

tssampid <- scan(tsidfile,"character")
tsgeneid <- scan(tsgenefile, "character")
tsexp <- scan(tsexpfile)
tsexp <- matrix(tsexp, ncol = length(tsgeneid), byrow=TRUE)
rownames(tsexp) <- tssampid
colnames(tsexp) <- tsgeneid

#make individual TS exp files for the 9 tissues.
tsexp <- cbind(info[,1:3],tsexp)
tislist <- scan(my.dir %&% "nine.spaces.tissue.list", "c", sep='\n')
for(tis in tislist){
  tisexp <- filter(tsexp,TISSUE==tis)
  write.table(tisexp,file=otd.dir %&% "obs_otd_exp/TS_" %&% tis %&% "_exp.txt",quote=F,row.names=F,sep="\t")
}
