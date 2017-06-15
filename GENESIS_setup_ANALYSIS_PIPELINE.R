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
mypcrel <- pcrelate(genoData = genoData, pcMat = mypcair$vectors[,1:20],training.set = mypcair$unrels)

pheno <- as.vector(as.matrix(read.table(phenoFile,header=F,na.string="NA")['V2']))
pheno <- pheno - 1

scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID = mypcrel$sample.id,pc1 = mypcair$vectors[,1],pc2 = mypcair$vectors[,2],
	pc3 = mypcair$vectors[,3],pc4 = mypcair$vectors[,4],pc5 = mypcair$vectors[,5],pc6 = mypcair$vectors[,6],
	pc7 = mypcair$vectors[,7],pc8 = mypcair$vectors[,8], pc9 = mypcair$vectors[,9],
	pc10 = mypcair$vectors[,10],pc11 = mypcair$vectors[,11],pc12 = mypcair$vectors[,12],
	pc13 = mypcair$vectors[,13],pc14 = mypcair$vectors[,14],pc15 = mypcair$vectors[,15],
	pc16 = mypcair$vectors[,16],pc17 = mypcair$vectors[,17],pc18 = mypcair$vectors[,18],
	pc19 = mypcair$vectors[,19],pc20 = mypcair$vectors[,20],pheno = pheno))
covMatList <- list("Kin" = pcrelateMakeGRM(mypcrel))

# creates a binary file -- can open in R with load()
save(scanAnnot, covMatList, mypcair, file = paste(filePrefix, '_GENESIS', sep=""))