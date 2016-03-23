####by Heather E. Wheeler 20141113####
"%&%" = function(a,b) paste(a,b,sep="")

###protein coding genes were extracted from the GTEx RNA-Seq data
###on genegate: /nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/2_pull_protein_coding_RPKM_rebuild.pl
###and scp'd to tarbell: /group/im-lab/hwheeler/cross-tissue/gtex-rnaseq/

rna.dir <- "/group/im-lab/hwheeler/cross-tissue/gtex-rnaseq/"
annot.dir <- "/group/im-lab/hwheeler/cross-tissue/gtex-annot/"
gt.dir <- "/group/im-lab/hwheeler/cross-tissue/gtex-genotypes/"

rnaid <- read.table(rna.dir %&% 'GTEx_Analysis_2014-06-13.RNA-seq.ID.list')
subid <- read.table(annot.dir %&% 'GTEx_Data_2014-06-13_Annotations_SubjectSampleMappingDS.txt',header=T)
tisid <- read.table(annot.dir %&% 'GTEx_Analysis_2014-06-13.SampleTissue.annot',header=T,sep="\t")
sexid <- read.table(annot.dir %&% 'GTEx_Data_2014-06-13_Annotations_SubjectID.GENDER.txt',header=T,sep="\t")
pcsid <- read.table(gt.dir %&% 'GTEx_Analysis_2014-06-13_OMNI_2.5M_5M_451Indiv_PostImput_20genotPCs.txt',header=T)
pcs <- pcsid[,1:5]

a<-merge(rnaid,subid,by.x="V1",by.y="SAMPID",sort=F)
b<-merge(a,tisid,by.x="V1",by.y="SAMPID",sort=F)
c<-merge(b,sexid,by.x="SUBJID",by.y="SUBJID",sort=F)
d<-merge(c,pcs,by.x="SUBJID",by.y="FID",sort=F,all.x=T)
e<-merge(rnaid,d,by.x="V1",by.y="V1",sort=F)
f<-e[,-6]
final <- f[,-3]

write.table(final,file=rna.dir %&% "GTEx_Analysis_2014-06-13.RNA-seq.IDinfo.list",row.names=F,col.names=F,quote=F,sep="\t")

##SAMPID SUBJID SMTSD GENDER PC1 PC2 PC3##
