###get obs v sim correlations for plotting
  
library(dplyr)
library(data.table)
library(tidyr)
"%&%" = function(a,b) paste(a,b,sep="")
pre.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
my.dir <- pre.dir %&% "paper-reviewer-requests/OTD_simulations/"
obs.dir <- my.dir %&% "obs_otd_exp/"
sim.dir <- my.dir %&% "lmer.fits/"
obs.h2.dir <- pre.dir %&% "gtex-h2-estimates/"

###read in obs OTD
obsct <- fread(obs.dir %&% "CT_exp.txt")
tislist <- scan(pre.dir %&% "nine.spaces.tissue.list", "c", sep='\n')
for(tis in tislist){
  obstis <- fread(obs.dir %&% "TS_" %&% tis %&% "_exp.txt", sep='\t')
  if(exists("obsts") == FALSE){
    obsts <- obstis
  }else{
    obsts <- rbind(obsts, obstis)
  }
}

###compare CT sim vs obs
simlist <- scan(my.dir %&% "shortsimlist", "c")
simfilelist <- scan(my.dir %&% "ctsimfilelist", "c")
for(i in 1:length(simlist)){
  simct <- fread(my.dir %&% simfilelist[i])
  s <- simlist[i]
  obsctdf <- data.frame(obsct) %>% filter(SUBJID %in% simct$SUBJID) %>% select(-SUBJID)
  simctdf <- data.frame(simct) %>% select(-SUBJID)
  corbygene <- cor(obsctdf,simctdf)
  genecor <- diag(corbygene)
  genenames <- names(genecor)
  cordf <- mutate(data.frame(genecor),sim=s,gene=genenames)
  if(exists("allctcor") == FALSE){
    allctcor <- cordf
  }else{
    allctcor <- rbind(allctcor,cordf)
  }
}

write.table(allctcor, file=my.dir %&% "CT_genecor_obs_v_sim.txt",quote=F,row.names = F)
rm('allctcor')

###compare TS sim vs obs
tisannot <- read.table(pre.dir %&% "gtex-annot/GTEx_Analysis_2014-06-13.SampleTissue.annot",header=T,sep='\t')
simlist <- scan(my.dir %&% "shortsimlist", "c")
simfilelist <- scan(my.dir %&% "tssimfilelist", "c")
for(i in 1:length(simlist)){
  simts <- fread(sim.dir %&% simfilelist[i])
  s <- simlist[i]
  obstsdf <- data.frame(obsts) %>% filter(SAMPID %in% simts$SAMPID)
  simtsdf <- data.frame(simts)
  simtstis <- left_join(simtsdf,tisannot,by='SAMPID')
  for(tis in tislist){
    tissim <- dplyr::filter(simtstis,SMTSD==tis)
    tisobs <- dplyr::filter(obstsdf,TISSUE==tis)
    sim <- dplyr::select(tissim,-SAMPID,-SMTSD)
    obs <- dplyr::select(tisobs,-SAMPID,-SUBJID,-TISSUE)
    tscorbygene <- cor(obs,sim)
    tsgenecor <- diag(tscorbygene)
    genenames <- names(tsgenecor)
    cordf <- mutate(data.frame(tsgenecor),tissue=tis,sim=s,gene=genenames)
    write.table(cordf, file=my.dir %&% "working_TS_genecor_obs_v_sim.txt",append=T,quote=F,row.names=F,col.names=F,sep="\t")
    if(exists("allcor") == FALSE){
      allcor <- cordf
    }else{
      allcor <- rbind(allcor,cordf)
    }
  }
}

write.table(allcor, file=my.dir %&% "TS_genecor_obs_v_sim.txt",quote=F,row.names = F,sep="\t")
