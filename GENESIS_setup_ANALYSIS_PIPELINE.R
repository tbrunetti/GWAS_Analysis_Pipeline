#source("https://bioconductor.org/biocLite.R")
#biocLite("GWASTools")

# read in command line arugments -- POSITIONAL!
args <- commandArgs(trailingOnly = T)
filePrefix <- args[1] # for KING and PLINK files
phenoFile <- args[2] # specially designed for GENESIS (essentially the 2nd and 6th column of fam file; need full name and path)

library(GWASTools, lib.loc=args[3])

###Read in Barbados genotype data for PCAs
#biocLite("SNPRelate")
library("SNPRelate", lib.loc=args[3])
snpgdsBED2GDS(bed.fn = paste(filePrefix, ".bed", sep=""), bim.fn = paste(filePrefix, ".bim", sep=""), fam.fn = paste(filePrefix, ".fam", sep=""), 
              out.gdsfn = paste(filePrefix, ".gds", sep=""))

#biocLite("GENESIS")
library("GENESIS", lib.loc=args[3])
file.kin0 <- paste(filePrefix, ".kin0", sep="")
file.kin <- paste(filePrefix, ".kin", sep="")
geno <- GdsGenotypeReader(filename = paste(filePrefix, ".gds", sep=""))
genoData <- GenotypeData(geno)
iids <- getScanID(genoData)

###Run PC analysis
Kingmat <- king2mat(file.kin0=file.kin0,file.kin=NULL,type="kinship",iids = iids)
mypcair <- pcair(genoData = genoData, kinMat = Kingmat,divMat = Kingmat)
mypcrel <- pcrelate(genoData = genoData, pcMat = mypcair$vectors[,1:8],training.set = mypcair$unrels)

pheno <- as.vector(as.matrix(read.table(phenoFile,header=F,na.string="NA")['V2']))
pheno <- pheno - 1

scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID = mypcrel$sample.id,pc1 = mypcair$vectors[,1],pc2 = mypcair$vectors[,2],
	pc3 = mypcair$vectors[,3],pc4 = mypcair$vectors[,4],pc5 = mypcair$vectors[,5],pc6 = mypcair$vectors[,6],
	pc7 = mypcair$vectors[,7],pc8 = mypcair$vectors[,8], pheno = pheno))
covMatList <- list("Kin" = pcrelateMakeGRM(mypcrel))

# creates a binary file -- can open in R with load()
save(scanAnnot, covMatList, mypcair, file = paste(filePrefix, '_GENESIS', sep=""))