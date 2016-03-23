library(dplyr)
"%&%" = function(a,b) paste(a,b,sep="")

##combine tissue-specific pve50 results from BSLMM

my.dir <- '/group/im-lab/nas40t2/hwheeler/cross-tissue/'
bslmm.dir <- my.dir %&% 'BSLMM_exp/'
tislist <- scan(my.dir %&% 'GTEx_nine_tissues','c')

ct <- read.table(bslmm.dir %&% 'cross-tissue_exp_BSLMM-s100K_iterations_all_chr1-22_2015-08-06.txt',header=T)
allpve <- ct %>% select(gene,pve50)
colnames(allpve) <- c('gene','CrossTissue')

for(i in 1:length(tislist)){
  tis <- tislist[i]
  tisdata <- read.table(bslmm.dir %&% tis %&% '_TS_exp_BSLMM-s100K_iterations_all_chr1-22_2015-08-06.txt',header=T)
  tispve <- tisdata %>% select(gene,pve50)
  allpve <- left_join(allpve,tispve,by='gene')
  tis <- gsub('(','',tis,fixed=T) #replace any parentheses in tis string
  tis <- gsub(')','',tis,fixed=T)
  tis <- gsub('-','',tis)
  colnames(allpve)[2+i] <- tis
}

write.table(allpve,file="GTEx_Tissue-Specific_local_PVE_by_BSLMM.txt",quote=F,row.names=F)

##combine tissue-Wide pve50 results from BSLMM

my.dir <- '/group/im-lab/nas40t2/hwheeler/cross-tissue/'
bslmm.dir <- my.dir %&% 'BSLMM_exp/'
tislist <- scan(my.dir %&% 'GTEx_nine_tissues','c')

ct <- read.table(bslmm.dir %&% 'cross-tissue_exp_BSLMM-s100K_iterations_all_chr1-22_2015-08-06.txt',header=T)
allpve <- ct %>% select(gene,pve50)
colnames(allpve) <- c('gene','CrossTissue')

for(i in 1:length(tislist)){
  tis <- tislist[i]
  tisdata <- read.table(bslmm.dir %&% tis %&% '_TW_exp_BSLMM-s100K_iterations_all_chr1-22_2015-10-18.txt',header=T)
  tispve <- tisdata %>% select(gene,pve50)
  allpve <- left_join(allpve,tispve,by='gene')
  tis <- gsub('(','',tis,fixed=T) #replace any parentheses in tis string
  tis <- gsub(')','',tis,fixed=T)
  tis <- gsub('-','',tis)
  colnames(allpve)[2+i] <- tis
}

write.table(allpve,file="GTEx_Tissue-Wide_local_PVE_by_BSLMM.txt",quote=F,row.names=F)
