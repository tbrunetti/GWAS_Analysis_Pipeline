#! Rscript --vanilla --default-packages=utils
library("gdsfmt", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3/")
library("SNPRelate", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3/")
library("GWASTools", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3/")
library("SeqVarTools", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3/") 
library("GENESIS", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3/")

args <- commandArgs(trailingOnly = T)
chunkedFile = args[1]
pathPLINKprefix = args[2]

load(paste(pathPLINKprefix, "_GENESIS", sep=""))
###association testing
snpgdsBED2GDS(bed.fn = paste(chunkedFile, ".bed", sep=""), bim.fn = paste(chunkedFile, '.bim', sep=""), fam.fn = paste(chunkedFile, ".fam", sep=""), out.gdsfn = paste(chunkedFile, ".gds", sep=""))

gds <- GdsGenotypeReader(paste(chunkedFile, ".gds", sep=""))
genoData <- GenotypeData(gds)

nullmod.bin <- fitNullMM(scanData = scanAnnot, outcome = "pheno", covars = c("pc1", "pc8"), covMatList = covMatList, family=binomial(link = "logit"))
myassoc <- assocTestMM(genoData = genoData, nullMMobj = nullmod.bin, test="Score")
write.table(myassoc,file=paste(chunkedFile,".results.txt",sep=''),sep="\t",row.names=F,col.names=T,quote=F)
