####by Heather E. Wheeler 20160617####
"%&%" = function(a,b) paste(a,b,sep="")
args <- commandArgs(trailingOnly=T)
date = Sys.Date()
library(dplyr)
library(data.table)
my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
lmer.dir <- my.dir %&% "lmer.fits/"
annot.dir <- my.dir %&% "gtex-annot/"
otd.dir <- my.dir %&% "paper-reviewer-requests/OTD_simulations/"
sim.dir <- otd.dir %&% "sim_exp/"

###############################################
errvar = args[1] ## choose 'ct' or 'ts' or 'sum'
mult = as.numeric(args[2]) ##constant to multiply the variance by
seed = 123

errfunc <- function(x){
  n=length(x)
  v=var(x)*mult
  err=rnorm(n,mean=0,sd=v)
  return(err)
}

set.seed(seed)

ct <- fread(otd.dir %&% "obs_otd_exp/CT_exp.txt")
tislist <- scan(my.dir %&% "nine.spaces.tissue.list", "c", sep='\n')

for(tis in tislist){
  tisexp <- fread(otd.dir %&% "obs_otd_exp/TS_" %&% tis %&% "_exp.txt",sep='\t')

  ctexp <- data.frame(filter(ct, SUBJID %in% tisexp$SUBJID))
  tsexp <- data.frame(filter(tisexp, SUBJID %in% ctexp$SUBJID))

  ctmat <- as.matrix(ctexp[,2:dim(ctexp)[2]])
  tsmat <- as.matrix(tsexp[,4:dim(tsexp)[2]])

  presim <- ctmat + tsmat

  if(errvar == 'ct'){
    errmat <- apply(ctmat,2,errfunc) #generate error matrix from variance of each gene's CT exp
  }else if(errvar =='ts'){
    errmat <- apply(tsmat,2,errfunc) #generate error matrix from variance of each gene's TS exp
  }else{
    errmat <- apply(presim,2,errfunc) #generate error matrix from variance of each gene's CT+TS exp
  }

  sim <- ctmat + tsmat + errmat
  simout <- data.frame(sim)
  simout <- cbind(tsexp[,1:2], simout)

  if(exists("simall") == FALSE){
    simall <- simout
  }else{
    simall <- rbind(simall, simout)
  }
}

write.table(simall, sim.dir %&% "sim_exp_phenotype_errvar-" %&% errvar %&% "_mult-" %&% mult %&% 
              "_seed-" %&% seed %&% ".txt",quote=F, row.names=F)
