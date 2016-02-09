####by Heather E. Wheeler 20141120####
args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()
library(dplyr)

thresh <- 'fdr0.05' ##threshold for inclusion of trans FHS SNPs
tis <- "DGN-WB"

pre.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
my.dir <- pre.dir %&% "expArch_DGN-WB_imputedGTs/"
annot.dir <- pre.dir %&% "gtex-annot/"
grm.dir <- my.dir %&% "dgn-grms/FHS_eQTLs_" %&% thresh %&% "/"
exp.dir <- my.dir %&% "dgn-exp/"

gencodefile <- annot.dir %&% "gencode.v18.genes.patched_contigs.summary.protein.chr" %&% args[1]
gencodeset <- args[1]

###############################################
### Scan expression data
idfile <- exp.dir %&% "DGN-WBexp.ID.list"
genefile <- exp.dir %&% "DGN-WBexp.GENE.list"
expfile <- exp.dir %&% "DGN-WB.rntransform.exp.IDxGENE"

subjid <- scan(idfile,"character")
geneid <- scan(genefile, "character")
expdata <- scan(expfile)
expdata <- matrix(expdata, ncol = length(geneid), byrow=TRUE)
rownames(expdata) <- subjid
colnames(expdata) <- geneid

### Get gene subset to analyze
system("ls " %&% grm.dir %&% "local*id > tmp." %&% thresh %&% gencodeset)
system("awk -F \"local-\" \'{print $2}\' < tmp." %&% thresh %&% gencodeset %&% " | awk -F \".grm\" \'{print $1}\' >localGRM.list." %&% thresh %&% gencodeset)
localfile <- "localGRM.list."	%&% thresh %&% gencodeset
locallist <- as.data.frame(scan(localfile,"character")) ##list of local GRMs
colnames(locallist)<-'gene'

gencode <- read.table(gencodefile)
rownames(gencode) <- gencode[,5]
colnames(gencode) <- c('chr','str','start','end','ensid','gene','func','known')

grmlist <- inner_join(locallist,gencode,by='gene') ##genes in genecodeset with grm
geneidlist <- as.data.frame(geneid)
colnames(geneidlist) <- 'gene'
ensidlist <- inner_join(grmlist,geneidlist,by='gene') ##genes with grm and exp
rownames(ensidlist) <- ensidlist$ensid
ensgene <- as.character(ensidlist$gene)

### Get individual subset to analyze
grm.id <- read.table(grm.dir %&% "DGN.global_Chr" %&% gencodeset %&% ".grm.id")
indidlist <- intersect(subjid,grm.id[,1]) ##subjects with exp and grm
nsubj <- length(indidlist)

### Get expression data from intersected genes and subjects
localexp <- expdata[indidlist,ensgene]
localensid <- ensidlist$ensid

### Output matrices
loc.mat <- matrix(0,nrow=length(localensid),ncol=7)
colnames(loc.mat) <- c("tissue","N","ensid","gene","local.h2","local.se","local.p")

glo.mat <- matrix(0,nrow=length(localensid),ncol=4)
colnames(glo.mat) <- c("gene","global.h2","global.se","global.p")

locglo.mat <- matrix(0,nrow=length(localensid),ncol=5)
colnames(locglo.mat) <- c("gene","loc.jt.h2","loc.jt.se","glo.jt.h2","glo.jt.se")

for(i in 1:length(localensid)){
	cat(i,"of",length(localensid),"\n")
	ensid <- as.character(localensid[i])
	gene <- as.character(gencode[ensid,6])
	chr <- as.character(gencode[ensid,1])
	c <- as.numeric(substr(chr,4,5))

	#output expression pheno for gcta
	geneexp <- cbind(rownames(localexp),localexp[,gene])
	write.table(geneexp, file="tmp.pheno." %&% thresh %&% gencodeset, col.names=F, quote=F) #output pheno for gcta

	## Y ~ localGRM
	runLOC <- "gcta64 --grm " %&% grm.dir %&% "local-" %&% gene %&% " --reml-no-constrain --pheno tmp.pheno." %&% thresh %&% gencodeset %&% " --out tmp." %&% thresh %&% gencodeset
	system(runLOC)
	if(file.exists("tmp." %&% thresh %&% gencodeset %&% ".hsq")==TRUE){
		hsq <- scan("tmp." %&% thresh %&% gencodeset %&% ".hsq","character")
        	res <- c(tis, nsubj, ensid, gene, hsq[14], hsq[15], hsq[25])
	}else{ #gcta did not coverge or Error: the information matrix is not invertible.
		res <- c(tis, nsubj, ensid, gene, NA, NA, NA)
	}
	loc.mat[i,] <- res
	system("rm tmp." %&% thresh %&% gencodeset %&% ".hsq")

	## Y ~ globalGRM (all FHS eQTLs, not including localGRM SNPs) ##changed 3/18/15 to all FHS eQTLs on other chrs
	###combine GRMs
	chrlist <- c(1:22)
	otherchrs <- setdiff(chrlist,c)
	grmmat <- matrix(NA,nrow=length(otherchrs),ncol=1)
	for(j in 1:length(otherchrs)){
		chrom <- otherchrs[j]
		grmfile <- grm.dir %&% "DGN.global_Chr" %&% chrom
		grmmat[j,] <- grmfile
	}
	write.table(grmmat,file="tmp.grm.list" %&% thresh %&% gencodeset,quote=F,col.names=F,row.names=F)
	runGRM <- "gcta64 --mgrm-bin tmp.grm.list" %&% thresh %&% gencodeset %&% " --make-grm --out " %&% grm.dir %&% "global-" %&% gene %&% "-otherChr"
	system(runGRM)

	runGLO <- "gcta64 --grm " %&% grm.dir %&% "global-" %&% gene %&% "-otherChr --reml-no-constrain --pheno tmp.pheno." %&% thresh %&% gencodeset %&% " --out tmp." %&% thresh %&% gencodeset
	system(runGLO)
	if(file.exists("tmp." %&% thresh %&% gencodeset %&% ".hsq")==TRUE){
		hsq <- scan("tmp." %&% thresh %&% gencodeset %&% ".hsq","character")
		res <- c(gene, hsq[14], hsq[15], hsq[25])
	}else{
                res <- c(gene, NA, NA, NA) ##gcta did not converge
        }
	glo.mat[i,] <- res 
	system("rm tmp." %&% thresh %&% gencodeset %&% ".hsq")

	## Y ~ localGRM + globalGRM (joint model)
	runMULT <- "echo " %&% grm.dir %&% "local-" %&% gene %&% " > tmp.multiGRM." %&% thresh %&% gencodeset
	runMULT2 <- "echo " %&% grm.dir %&% "global-" %&% gene %&% "-otherChr >> tmp.multiGRM." %&% thresh %&% gencodeset
	system(runMULT)
        system(runMULT2)
	runLOCGLO <- "gcta64 --mgrm-bin tmp.multiGRM." %&% thresh %&% gencodeset %&% " --reml-no-constrain --pheno tmp.pheno." %&% thresh %&% gencodeset %&% " --out tmp." %&% thresh %&% gencodeset
	system(runLOCGLO)
	if(file.exists("tmp." %&% thresh %&% gencodeset %&% ".hsq")==TRUE){
		hsq <- scan("tmp." %&% thresh %&% gencodeset %&% ".hsq","character")
		res <- c(gene, hsq[17], hsq[18], hsq[20], hsq[21])
	}else{
		res <- c(gene, NA, NA, NA, NA) ##gcta did not converge
	}
	locglo.mat[i,] <- res
	system("rm tmp." %&% thresh %&% gencodeset %&% ".hsq")
}


full.mat <- cbind(loc.mat,glo.mat,locglo.mat)
write.table(full.mat,file=tis %&% ".h2.all.models_FHS" %&% thresh %&% ".Chr" %&% gencodeset %&% "_globaleQTLOtherChr_reml-no-constrain." %&% date %&% ".txt",quote=F,row.names=F,sep="\t")
