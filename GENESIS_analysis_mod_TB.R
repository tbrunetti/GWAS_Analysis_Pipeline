args <- commandArgs(trailingOnly = T)
chunkedFile = args[1]
pathPLINKprefix = args[2]
rlib = args[3]

#! Rscript --vanilla --default-packages=utils
library("gdsfmt", lib.loc=rlib)
library("SNPRelate", lib.loc=rlib)
library("GWASTools", lib.loc=rlib)
library("SeqVarTools", lib.loc=rlib) 
library("GENESIS", lib.loc=rlib)


load(paste(pathPLINKprefix, "_GENESIS", sep=""))
###association testing
snpgdsBED2GDS(bed.fn = paste(chunkedFile, ".bed", sep=""), bim.fn = paste(chunkedFile, '.bim', sep=""), fam.fn = paste(chunkedFile, ".fam", sep=""), out.gdsfn = paste(chunkedFile, ".gds", sep=""))

gds <- GdsGenotypeReader(paste(chunkedFile, ".gds", sep=""))
genoData <- GenotypeData(gds)

nullmod.bin <- fitNullMM(scanData = scanAnnot, outcome = "pheno", covars = c("pc1", "pc8"), covMatList = covMatList, family=binomial(link = "logit"))
myassoc <- assocTestMM(genoData = genoData, nullMMobj = nullmod.bin, test="Score")
write.table(myassoc,file=paste(chunkedFile,".results.txt",sep=''),sep="\t",row.names=F,col.names=T,quote=F)
