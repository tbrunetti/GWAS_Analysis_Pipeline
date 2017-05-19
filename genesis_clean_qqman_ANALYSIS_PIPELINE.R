library("qqman", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3/")

args <- commandArgs(trailingOnly = T)

final_merged_text = args[1]
final_bim = args[2]

dat.info<-read.delim(final_merged_text)
bim_file_for_merged_data <- read.delim(final_bim, col.names = c('chr', 'snpID', 'pos_centMorgans', 'POS', 'allele_1', 'allele_2'))

colnames(dat.info)
dim(dat.info)


dat.info.out <- merge(dat.info, bim_file_for_merged_data, by=c('snpID', 'chr'), all.y = T)



#common SNPs
dat.info.out.common<-dat.info.out[which(dat.info.out$MAF >= 0.05),]
dim(dat.info.out.common)
summary(dat.info.out.common$Score.pval)

#rare SNPs
dat.info.out.rare<-dat.info.out[which(dat.info.out$MAF < 0.05),]
dim(dat.info.out.rare)
summary(dat.info.out.rare$Score.pval)

#all qq plot
observed <- sort(dat.info.out$Score.pval)
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
bmp(filename="qqplot_typed_overlap_allcohort_new.bmp", width=800, height=800, bg="white", type="cairo")
#pdf("qqplot_typed_overlap_allcohort.pdf", width=6, height=6)
plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
dev.off()


#inflation factor lambda
chisq2 <- qchisq(1-dat.info.out$Score.pval,1,lower.tail = T)
median(chisq2,na.rm=T)/qchisq(0.5,1)#lambda1

#manhattan plot
library(qqman)
bmp(filename="manhattan_typed_overlap_allcohort_new.bmp", width=800, height=600, bg="white", type="cairo")
#pdf("manhattan_typed_overlap_allcohort.pdf",width=21,height=10)
par(font.axis = 2)
manhattan(dat.info.out,chr = "chr", bp = "POS", p = "Score.pval", snp = "snpID",col = c("gray60", "gray10"), chrlabs = NULL,highlight = NULL, logp = TRUE,suggestiveline = F, genomewideline = F,ylim=c(0,10),main="Association analysis:Genesis Washington")
dev.off()

#common qq plot
observed <- sort(dat.info.out.common$Score.pval)
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
bmp(filename="qqplot_typed_overlap_allcohort_common_new.bmp", width=800, height=800, bg="white", type="cairo")
#pdf("qqplot_typed_overlap_allcohort_common.pdf", width=6, height=6)
plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
dev.off()

#inflation factor lambda
chisq2 <- qchisq(1-dat.info.out.common$Score.pval,1,lower.tail = T)
median(chisq2,na.rm=T)/qchisq(0.5,1)#lambda1

#manhattan plot common
library(qqman)
bmp(filename="manhattan_typed_overlap_allcohort_common_new.bmp", width=800, height=600, bg="white", type="cairo")
#pdf("manhattan_typed_overlap_allcohort_common_new.pdf",width=21,height=10)
par(font.axis = 2)
manhattan(dat.info.out.common,chr = "chr", bp = "POS", p = "Score.pval", snp = "snpID",col = c("gray60", "gray10"), chrlabs = NULL,highlight = NULL, logp = TRUE,suggestiveline = F, genomewideline = F, ylim=c(0,10) ,main="Association analysis:Genesis Washington common variants")
dev.off()

#rare
observed <- sort(dat.info.out.rare$Score.pval)
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
bmp(filename="qqplot_typed_overlap_allcohort_rare_new.bmp", width=800, height=800, bg="white", type="cairo")
#pdf("qqplot_typed_overlap_allcohort_rare.pdf", width=6, height=6)
plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
dev.off()

#inflation factor lambda
chisq2 <- qchisq(1-dat.info.out.rare$Score.pval,1,lower.tail = T)
median(chisq2,na.rm=T)/qchisq(0.5,1)#lambda1

#manhattan plot rare
library(qqman)
bmp(filename="manhattan_typed_overlap_allcohort_rare_new.bmp", width=800, height=600, bg="white", type="cairo")
#pdf("manhattan_typed_overlap_allcohort_rare.pdf",width=21,height=10)
par(font.axis = 2)
manhattan(dat.info.out.rare,chr = "chr", bp = "POS", p = "Score.pval", snp = "snpID",col = c("gray60", "gray10"), chrlabs = NULL,highlight = NULL, logp = TRUE,suggestiveline = F, genomewideline = F, ylim=c(0,10) ,main="Association analysis:Genesis Washington rare variants")
#dev.off()


#qq plot all snps with CI

bmp(filename="qqplot_typed_overlap_allcohort_new_with_CI.bmp", width=800, height=800, bg="white", type="cairo")

## obs <- readfile; p-values only
## read in your p-values,
## here I generated some
obs<- dat.info.out$Score.pval
N <- 1000000 ## number of p-values
## create the null distribution
## (-log10 of the uniform)
null <- -log(1:N/N,10)
MAX <- max(c(obs,null))
## create the confidence intervals
c95 <- rep(0,N)
c05 <- rep(0,N)
## the jth order statistic from a
## uniform(0,1) sample
## has a beta(j,n-j+1) distribution
## (Casella & Berger, 2002,
## 2nd edition, pg 230, Duxbury)
for(i in 1:N){
  c95[i] <- qbeta(0.95,i,N-i+1)
  c05[i] <- qbeta(0.05,i,N-i+1)
}
## plot the two confidence lines
plot(null, -log(c95,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l",
     axes=FALSE, xlab="", ylab="")
par(new=T)
plot(null, -log(c05,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l",
     axes=FALSE, xlab="", ylab="")
## add the diagonal
abline(0,1,col="red")
par(new=T)

observed <- sort(dat.info.out$Score.pval)
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
#pdf("qqplot_typed_overlap_allcohort.pdf", width=6, height=6)
plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
dev.off()

chisq2 <- qchisq(1-dat.info.out$Score.pval,1,lower.tail = T)
median(chisq2,na.rm=T)/qchisq(0.5,1)#lambda1


# qq common with CI
bmp(filename="qqplot_typed_overlap_allcohort_common_new_with_CI.bmp", width=800, height=800, bg="white", type="cairo")

obs<- dat.info.out.common$Score.pval
N <- 1000000 ## number of p-values
## create the null distribution
## (-log10 of the uniform)
null <- -log(1:N/N,10)
MAX <- max(c(obs,null))
## create the confidence intervals
c95 <- rep(0,N)
c05 <- rep(0,N)
## the jth order statistic from a
## uniform(0,1) sample
## has a beta(j,n-j+1) distribution
## (Casella & Berger, 2002,
## 2nd edition, pg 230, Duxbury)
for(i in 1:N){
  c95[i] <- qbeta(0.95,i,N-i+1)
  c05[i] <- qbeta(0.05,i,N-i+1)
}
## plot the two confidence lines
plot(null, -log(c95,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l",
     axes=FALSE, xlab="", ylab="")
par(new=T)
plot(null, -log(c05,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l",
     axes=FALSE, xlab="", ylab="")
## add the diagonal
abline(0,1,col="red")
par(new=T)


observed <- sort(dat.info.out.common$Score.pval)
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
#pdf("qqplot_typed_overlap_allcohort_common.pdf", width=6, height=6)
plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
dev.off()

chisq2 <- qchisq(1-dat.info.out.common$Score.pval,1,lower.tail = T)
median(chisq2,na.rm=T)/qchisq(0.5,1)#lambda1

# qq plot rare with CI
bmp(filename="qqplot_typed_overlap_allcohort_rare_new_with_CI.bmp", width=800, height=800, bg="white", type="cairo")

obs<- dat.info.out.rare$Score.pval
N <- 1000000 ## number of p-values
## create the null distribution
## (-log10 of the uniform)
null <- -log(1:N/N,10)
MAX <- max(c(obs,null))
## create the confidence intervals
c95 <- rep(0,N)
c05 <- rep(0,N)
## the jth order statistic from a
## uniform(0,1) sample
## has a beta(j,n-j+1) distribution
## (Casella & Berger, 2002,
## 2nd edition, pg 230, Duxbury)
for(i in 1:N){
  c95[i] <- qbeta(0.95,i,N-i+1)
  c05[i] <- qbeta(0.05,i,N-i+1)
}
## plot the two confidence lines
plot(null, -log(c95,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l",
     axes=FALSE, xlab="", ylab="")
par(new=T)
plot(null, -log(c05,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l",
     axes=FALSE, xlab="", ylab="")
## add the diagonal
abline(0,1,col="red")
par(new=T)

observed <- sort(dat.info.out.rare$Score.pval)
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
#pdf("qqplot_typed_overlap_allcohort_rare.pdf", width=6, height=6)
plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
dev.off()

chisq2 <- qchisq(1-dat.info.out.rare$Score.pval,1,lower.tail = T)
median(chisq2,na.rm=T)/qchisq(0.5,1)#lambda1
