load('/home/tonya/Desktop/CAAPA_MEGA_passed_QC_final_converted_African_American_maf_greater_thresh_hetFiltered_dups_removed_thousGen_GENESIS')
get_dataframe <- pData(scanAnnot)
write.table(get_dataframe, file='/home/tonya/Desktop/get_dataframe.txt', sep = '\t')


pcs_with_pops <- read.table('/home/tonya/Desktop/merged_pop_BASS_CAAPA.txt')

#Plot first N Eigenvalues from mypcair
data <- data.frame(as.matrix(mypcair$values))
data['eigenvalues'] <- data.frame(as.matrix(mypcair$values))
data['proportion_var_explained'] <-data.frame((as.matrix(mypcair$values/mypcair$sum.values)))
data['PC'] <- data.frame((as.numeric(rownames(data))))
subset_data <- subset(data,PC<201)
subset_data['new_eigensum'] <-data.frame(colSums(subset_data)['eigenvalues'])
subset_data['subset_Prop_variance'] <- data.frame(subset_data$eigenvalues/subset_data$new_eigensum)
plot(subset_data$PC, subset_data$subset_Prop_variance, col="black", pch=20, cex=1,ylab="Proportion of Variance Explained",xlab="",xaxt="n", main="Proportion of Variance explained by PCs (eigenvalue)")
axis(1,at=subset_data$PC,labels=subset_data$PC)

# Assign experimental vs TGP
pcs_with_pops$group[pcs_with_pops$population=="BASS"] <- "BASS"
pcs_with_pops$group[pcs_with_pops$population!="BASS"] <- "TGP"
pcs_with_pops$group <- as.factor(pcs_with_pops$group)

# scatter PC1 vs PC2
cols=c("blue", "darkgreen", "orange", "magenta", "pink", "red")
plot(pcs_with_pops$pc1, pcs_with_pops$pc2, col=cols[as.numeric(pcs_with_pops$population)],
     pch=c(17,20)[as.numeric(pcs_with_pops$group)], cex=0.5,xlab="Principal Component 1",ylab="Principal Component 2")
legend("topright", c("AFR","AMR","BASS","EAS","EUR","SAS"),pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=1,cex=.5)
legend("bottomleft", c("BASS","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)

# scatter PC1 vs PC3
cols=c("blue", "darkgreen", "orange", "magenta", "pink", "red")
plot(pcs_with_pops$pc1, pcs_with_pops$pc3, col=cols[as.numeric(pcs_with_pops$population)],
     pch=c(17,20)[as.numeric(pcs_with_pops$group)], cex=0.5,xlab="Principal Component 1",ylab="Principal Component 3")
legend("topright", c("AFR","AMR","BASS","EAS","EUR","SAS"),pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=1,cex=0.5)
legend("bottomleft", c("BASS","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)

# scatter PC1 vs PC4
cols=c("blue", "darkgreen", "orange", "magenta", "pink", "red")
plot(pcs_with_pops$pc1, pcs_with_pops$pc4, col=cols[as.numeric(pcs_with_pops$population)],
     pch=c(17,20)[as.numeric(pcs_with_pops$group)], cex=0.5,xlab="Principal Component 1",ylab="Principal Component 4")
legend("topright", c("AFR","AMR","BASS","EAS","EUR","SAS"),pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=1,cex=0.5)
legend("bottomleft", c("BASS","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)


# scatter PC1 vs PC5
cols=c("blue", "darkgreen", "orange", "magenta", "pink", "red")
plot(pcs_with_pops$pc1, pcs_with_pops$pc5, col=cols[as.numeric(pcs_with_pops$population)],
     pch=c(17,20)[as.numeric(pcs_with_pops$group)], cex=0.5,xlab="Principal Component 1",ylab="Principal Component 5")
legend("topright", c("AFR","AMR","BASS","EAS","EUR","SAS"),pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=1,cex=0.5)
legend("bottomleft", c("BASS","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)


# scatter PC2 vs PC3
cols=c("blue", "darkgreen", "orange", "magenta", "pink", "red")
plot(pcs_with_pops$pc2, pcs_with_pops$pc3, col=cols[as.numeric(pcs_with_pops$population)],
     pch=c(17,20)[as.numeric(pcs_with_pops$group)], cex=0.5,xlab="Principal Component 2",ylab="Principal Component 3")
legend("topright", c("AFR","AMR","BASS","EAS","EUR","SAS"),pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=1,cex=0.5)
legend("bottomleft", c("BASS","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)

# scatter PC2 vs PC4
cols=c("blue", "darkgreen", "orange", "magenta", "pink", "red")
plot(pcs_with_pops$pc2, pcs_with_pops$pc4, col=cols[as.numeric(pcs_with_pops$population)],
     pch=c(17,20)[as.numeric(pcs_with_pops$group)], cex=0.5,xlab="Principal Component 2",ylab="Principal Component 4")
legend("topright", c("AFR","AMR","BASS","EAS","EUR","SAS"),pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=1,cex=0.5)
legend("bottomleft", c("BASS","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)

# scatter PC2 vs PC5
cols=c("blue", "darkgreen", "orange", "magenta", "pink", "red")
plot(pcs_with_pops$pc2, pcs_with_pops$pc5, col=cols[as.numeric(pcs_with_pops$population)],
     pch=c(17,20)[as.numeric(pcs_with_pops$group)], cex=0.5,xlab="Principal Component 2",ylab="Principal Component 5")
legend("topright", c("AFR","AMR","BASS","EAS","EUR","SAS"),pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=1,cex=0.5)
legend("bottomleft", c("BASS","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)

# scatter PC3 vs PC4
cols=c("blue", "darkgreen", "orange", "magenta", "pink", "red")
plot(pcs_with_pops$pc3, pcs_with_pops$pc4, col=cols[as.numeric(pcs_with_pops$population)],
     pch=c(17,20)[as.numeric(pcs_with_pops$group)], cex=0.5,xlab="Principal Component 3",ylab="Principal Component 4")
legend("topright", c("AFR","AMR","BASS","EAS","EUR","SAS"),pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=1,cex=0.5)
legend("bottomleft", c("BASS","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)


# scatter PC3 vs PC5
cols=c("blue", "darkgreen", "orange", "magenta", "pink", "red")
plot(pcs_with_pops$pc3, pcs_with_pops$pc5, col=cols[as.numeric(pcs_with_pops$population)],
     pch=c(17,20)[as.numeric(pcs_with_pops$group)], cex=0.5,xlab="Principal Component 3",ylab="Principal Component 5")
legend("topright", c("AFR","AMR","BASS","EAS","EUR","SAS"),pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=1,cex=0.5)
legend("bottomleft", c("BASS","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)

#Boxplots centered on mean African ancestry
PCA_pops_AFR = subset(pcs_with_pops,population=="AFR")
#boxplot for PC1
m11=mean(PCA_pops_AFR$pc1)
s11=sd(PCA_pops_AFR$pc1)
boxplot(pc11 ~ population, data = pcs_with_pops,ylab = "PC1", xlab="PC1, +/- 6 SD on AFR",cex=0.5,outline=F)
stripchart(pc1 ~ population, data = pcs_with_pops,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m11, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m11-(6*s11), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m11+(6*s11), col = "black", lty=2, lwd=1,xpd=F)

#boxplot for PC2
m12=mean(PCA_pops_AFR$pc2)
s12=sd(PCA_pops_AFR$pc2)
boxplot(pc2 ~ population, data = pcs_with_pops,ylab = "PC2", xlab="PC2, +/- 6 SD on AFR",cex=0.5,outline=F)
stripchart(pc2 ~ population, data = pcs_with_pops,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m12, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m12-(6*s12), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m12+(6*s12), col = "black", lty=2, lwd=1,xpd=F)

#boxplot for PC3
m13=mean(PCA_pops_AFR$pc3)
s13=sd(PCA_pops_AFR$pc3)
boxplot(pc3 ~ population, data = pcs_with_pops,ylab = "PC3", xlab="PC3, +/- 6 SD on AFR",cex=0.5,outline=F)
stripchart(pc3 ~ population, data = pcs_with_pops,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m13, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m13-(6*s13), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m13+(6*s13), col = "black", lty=2, lwd=1,xpd=F)

#boxplot for PC4
m14=mean(PCA_pops_AFR$pc4)
s14=sd(PCA_pops_AFR$pc4)
boxplot(pc4 ~ population, data = pcs_with_pops,ylab = "PC4", xlab="PC4, +/- 6 SD on AFR",cex=0.5,outline=F)
stripchart(pc4 ~ population, data = pcs_with_pops,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m14, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m14-(6*s14), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m14+(6*s14), col = "black", lty=2, lwd=1,xpd=F)

#boxplot for PC5
m15=mean(PCA_pops_AFR$pc5)
s15=sd(PCA_pops_AFR$pc5)
boxplot(pc5 ~ population, data = pcs_with_pops,ylab = "PC5", xlab="PC5, +/- 6 SD on AFR",cex=0.5,outline=F)
stripchart(pc5 ~ population, data = pcs_with_pops,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m15, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m15-(6*s15), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m15+(6*s15), col = "black", lty=2, lwd=1,xpd=F)


#Determine outliers
BASS_only = subset(pcs_with_pops,population=="BASS")
BASS_only$PC1_outlier[as.numeric(BASS_only$pc1)<(m11-(6*s11)) | as.numeric(BASS_only$pc1)>(m11+(6*s11))] <- 1
BASS_only$PC1_outlier[as.numeric(BASS_only$pc1)>=(m11-(6*s11)) & as.numeric(BASS_only$pc1)<=(m11+(6*s11))] <- 0
BASS_only$PC2_outlier[as.numeric(BASS_only$pc2)<(m12-(10*s12)) | as.numeric(BASS_only$pc2)>(m12+(10*s12))] <- 1
BASS_only$PC2_outlier[as.numeric(BASS_only$pc2)>=(m12-(10*s12)) & as.numeric(BASS_only$pc2)<=(m12+(10*s12))] <- 0
BASS_only$PC3_outlier[as.numeric(BASS_only$pc3)<(m13-(6*s13)) | as.numeric(BASS_only$pc3)>(m13+(6*s13))] <- 1
BASS_only$PC3_outlier[as.numeric(BASS_only$pc3)>=(m13-(6*s13)) & as.numeric(BASS_only$pc3)<=(m13+(6*s13))] <- 0
BASS_only$PC4_outlier[as.numeric(BASS_only$pc4)<(m14-(6*s14)) | as.numeric(BASS_only$pc4)>(m14+(6*s14))] <- 1
BASS_only$PC4_outlier[as.numeric(BASS_only$pc4)>=(m14-(6*s14)) & as.numeric(BASS_only$pc4)<=(m14+(6*s14))] <- 0
BASS_only$PC5_outlier[as.numeric(BASS_only$pc5)<(m15-(6*s15)) | as.numeric(BASS_only$pc5)>(m15+(6*s15))] <- 1
BASS_only$PC5_outlier[as.numeric(BASS_only$pc5)>=(m15-(6*s15)) & as.numeric(BASS_only$pc5)<=(m15+(6*s15))] <- 0

BASS_outliers = subset(BASS_only,PC1_outlier==1 | PC2_outlier==1 | PC3_outlier==1 | PC4_outlier==1 | PC5_outlier==1)
write.table(BASS_outliers,"BASS_outliers.txt",row.names=F,col.names=T,quote=F,sep="\t")

ADRN_outliers <- read.table("ADRN_outliers.txt",header=T,all=T,sep="\t")
ADRN_outliers2 = subset(ADRN_outliers,select=c("PATIENT","AFFSTAT","PC1_outlier","PC2_outlier"))

ADRN_outliers2$EXTRACOL <- 1
ADRN2 <- merge(ADRN_only,ADRN_outliers2,all=T)
ADRNfin <- ADRN2[is.na(ADRN2$EXTRACOL),]
ADRNfin$EXTRACOL <- NULL
ADRN_final = subset(ADRNfin,select=c("PATIENT","PC1","PC2","PC3","PC4","PC5"))
write.table(ADRN_final,"ADRN_WGS_779.txt",row.names=F,col.names=T,quote=F,sep="\t")


#-------------------------------------------------------------------------------
#Plot first N Eigenvalues
data=read.table("ADRN_783_TGP_reduced_eigenvalues.txt", header=F,na.string="NA")
names(data) <- "PC"
data$Factor <- paste0("PC",seq.int(nrow(data)))
data$Row <- seq.int(nrow(data))
data2 = subset(data,Row<11)
datax = subset(data2,select="PC")
data2$Prop_Var <- data2$PC/colSums(datax)

#Plot series of Eigenvectors
pdf("plot.PCA.ADRN.783.TGPreduced.pdf",onefile=T)
par(mfrow=c(1,1))
par(mar=c(5,4.5,0,2))
par(omi=c(.5,0,.5,0))
fig=c(3,9,1,4)/10

plot(data2$Row, data2$Prop_Var, col="black", pch=20, cex=1,ylab="Proportion of Variance Explained",xlab="",xaxt="n")
axis(1,at=data2$Row,labels=data2$Factor)

TGP_pops <- read.table("TGP_pops.txt", header=T,na.string="NA")
PCA <- read.table("ADRN_783_TGP_reduced_eigenvectors.txt",header=T,na.string="NA")

PCA_pops = merge(TGP_pops,PCA,by="PATIENT",all.y=T)
#PCA_pops2 = subset(PCA_pops,is.na(Popn))
#write.table(PCA_pops2,"ADRN_remove.792.txt",row.names=F,col.names=T,quote=F,sep="\t")
PCA_pops$group[PCA_pops$Popn=="ADRN"] <- "ADRN"
PCA_pops$group[PCA_pops$Popn!="ADRN"] <- "TGP"
PCA_pops$group <- as.factor(PCA_pops$group)
PCA_pops$check <- as.numeric(PCA_pops$Popn)

#Scatterplot PC1 vs PC2
par(mfrow=c(1,1))
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
par(omi=c(.5,0,.5,0))
fig=c(3,9,1,4)/10
cols=c("blue", "darkgreen", "orange", "magenta", "pink", "red")
plot(PCA_pops$PC1, PCA_pops$PC2, col=cols[as.numeric(PCA_pops$Popn)],
     pch=c(17,20)[as.numeric(PCA_pops$group)], cex=0.5,xlab="Principal Component 1",ylab="Principal Component 2")
legend("topright", inset=c(-0.2,0), legend=c("ADRN","TGP"),pch=c(17,20),pt.cex=1,cex=1)
legend("topright", inset=c(-0.2,0.15), legend=c("ADRN","AFR","AMR","EAS","EUR","SAS"),
       pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=1,cex=1)

#Scatterplot PC1 vs PC3
par(mfrow=c(2,2))
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
par(omi=c(.5,0,.5,0))
fig=c(3,9,1,4)/10
plot(PCA_pops$PC1, PCA_pops$PC3, col=cols[as.numeric(PCA_pops$Popn)],cex.axis=0.5,cex.lab=0.5,
     pch=c(17,20)[as.numeric(PCA_pops$group)], cex=0.5,xlab="Principal Component 1",ylab="Principal Component 3")
legend("topright", inset=c(-0.25,0), legend=c("ADRN","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)
legend("topright", inset=c(-0.25,0.15), legend=c("ADRN","AFR","AMR","EAS","EUR","SAS"),
       pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=0.5,cex=0.5)

#Scatterplot PC1 vs PC4
plot(PCA_pops$PC1, PCA_pops$PC4, col=cols[as.numeric(PCA_pops$Popn)],cex.axis=0.5,cex.lab=0.5,
     pch=c(17,20)[as.numeric(PCA_pops$group)], cex=0.5,xlab="Principal Component 1",ylab="Principal Component 4")
legend("topright", inset=c(-0.25,0), legend=c("ADRN","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)
legend("topright", inset=c(-0.25,0.15), legend=c("ADRN","AFR","AMR","EAS","EUR","SAS"),
       pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=0.5,cex=0.5)

#Scatterplot PC1 vs PC5
plot(PCA_pops$PC1, PCA_pops$PC5, col=cols[as.numeric(PCA_pops$Popn)],cex.axis=0.5,cex.lab=0.5,
     pch=c(17,20)[as.numeric(PCA_pops$group)], cex=0.5,xlab="Principal Component 1",ylab="Principal Component 5")
legend("topright", inset=c(-0.25,0), legend=c("ADRN","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)
legend("topright", inset=c(-0.25,0.15), legend=c("ADRN","AFR","AMR","EAS","EUR","SAS"),
       pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=0.5,cex=0.5)

#Scatterplot PC2 vs PC3
plot(PCA_pops$PC2, PCA_pops$PC3, col=cols[as.numeric(PCA_pops$Popn)],cex.axis=0.5,cex.lab=0.5,
     pch=c(17,20)[as.numeric(PCA_pops$group)], cex=0.5,xlab="Principal Component 2",ylab="Principal Component 3")
legend("topright", inset=c(-0.25,0), legend=c("ADRN","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)
legend("topright", inset=c(-0.25,0.15), legend=c("ADRN","AFR","AMR","EAS","EUR","SAS"),
       pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=0.5,cex=0.5)

#Scatterplot PC2 vs PC4
plot(PCA_pops$PC2, PCA_pops$PC4, col=cols[as.numeric(PCA_pops$Popn)],cex.axis=0.5,cex.lab=0.5,
     pch=c(17,20)[as.numeric(PCA_pops$group)], cex=0.5,xlab="Principal Component 2",ylab="Principal Component 4")
legend("topright", inset=c(-0.25,0), legend=c("ADRN","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)
legend("topright", inset=c(-0.25,0.15), legend=c("ADRN","AFR","AMR","EAS","EUR","SAS"),
       pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=0.5,cex=0.5)

#Scatterplot PC2 vs PC5
plot(PCA_pops$PC2, PCA_pops$PC5, col=cols[as.numeric(PCA_pops$Popn)],cex.axis=0.5,cex.lab=0.5,
     pch=c(17,20)[as.numeric(PCA_pops$group)], cex=0.5,xlab="Principal Component 2",ylab="Principal Component 5")
legend("topright", inset=c(-0.25,0), legend=c("ADRN","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)
legend("topright", inset=c(-0.25,0.15), legend=c("ADRN","AFR","AMR","EAS","EUR","SAS"),
       pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=0.5,cex=0.5)

#Scatterplot PC3 vs PC4
plot(PCA_pops$PC3, PCA_pops$PC4, col=cols[as.numeric(PCA_pops$Popn)],cex.axis=0.5,cex.lab=0.5,
     pch=c(17,20)[as.numeric(PCA_pops$group)], cex=0.5,xlab="Principal Component 3",ylab="Principal Component 4")
legend("topright", inset=c(-0.25,0), legend=c("ADRN","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)
legend("topright", inset=c(-0.25,0.15), legend=c("ADRN","AFR","AMR","EAS","EUR","SAS"),
       pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=0.5,cex=0.5)

#Scatterplot PC3 vs PC5
plot(PCA_pops$PC3, PCA_pops$PC5, col=cols[as.numeric(PCA_pops$Popn)],cex.axis=0.5,cex.lab=0.5,
     pch=c(17,20)[as.numeric(PCA_pops$group)], cex=0.5,xlab="Principal Component 3",ylab="Principal Component 5")
legend("topright", inset=c(-0.25,0), legend=c("ADRN","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)
legend("topright", inset=c(-0.25,0.15), legend=c("ADRN","AFR","AMR","EAS","EUR","SAS"),
       pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=0.5,cex=0.5)

#Boxplots centered on mean European ancestry
PCA_pops_EUR = subset(PCA_pops,Popn=="EUR")

par(mfrow=c(1,1))
par(mar=c(5,4.5,0,2))
par(omi=c(.5,0,.5,0))
fig=c(3,9,1,4)/10
m11=mean(PCA_pops_EUR$PC1)
s11=sd(PCA_pops_EUR$PC1)
boxplot(PC1 ~ Popn, data = PCA_pops,ylab = "PC1", xlab="PC1, +/- 6 SD on EUR",cex=0.5,outline=F)
stripchart(PC1 ~ Popn, data = PCA_pops,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m11, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m11-(6*s11), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m11+(6*s11), col = "black", lty=2, lwd=1,xpd=F)

m12=mean(PCA_pops_EUR$PC2)
s12=sd(PCA_pops_EUR$PC2)
boxplot(PC2 ~ Popn, data = PCA_pops,ylab = "PC2", xlab="PC2, +/- 6 SD on EUR",cex=0.5,outline=F)
stripchart(PC2 ~ Popn, data = PCA_pops,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m12, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m12-(6*s12), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m12+(6*s12), col = "black", lty=2, lwd=1,xpd=F)

m13=mean(PCA_pops_EUR$PC3)
s13=sd(PCA_pops_EUR$PC3)
boxplot(PC3 ~ Popn, data = PCA_pops,ylab = "PC3", xlab="PC3, +/- 6 SD on EUR",cex=0.5,outline=F)
stripchart(PC3 ~ Popn, data = PCA_pops,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m13, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m13-(6*s13), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m13+(6*s13), col = "black", lty=2, lwd=1,xpd=F)

m14=mean(PCA_pops_EUR$PC4)
s14=sd(PCA_pops_EUR$PC4)
boxplot(PC4 ~ Popn, data = PCA_pops,ylab = "PC4", xlab="PC4, +/- 6 SD on EUR",cex=0.5,outline=F)
stripchart(PC4 ~ Popn, data = PCA_pops,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m14, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m14-(6*s14), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m14+(6*s14), col = "black", lty=2, lwd=1,xpd=F)

m15=mean(PCA_pops_EUR$PC5)
s15=sd(PCA_pops_EUR$PC5)
boxplot(PC5 ~ Popn, data = PCA_pops,ylab = "PC5", xlab="PC5, +/- 6 SD on EUR",cex=0.5,outline=F,ylim=c(-0.06,0.08))
stripchart(PC5 ~ Popn, data = PCA_pops,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m15, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m15-(6*s15), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m15+(6*s15), col = "black", lty=2, lwd=1,xpd=F)

#Determine outliers
ADRN_only = subset(PCA_pops,Popn=="ADRN")
ADRN_only$PC1_outlier[as.numeric(ADRN_only$PC1)<(m11-(6*s11)) | as.numeric(ADRN_only$PC1)>(m11+(6*s11))] <- 1
ADRN_only$PC1_outlier[as.numeric(ADRN_only$PC1)>=(m11-(6*s11)) & as.numeric(ADRN_only$PC1)<=(m11+(6*s11))] <- 0
ADRN_only$PC2_outlier[as.numeric(ADRN_only$PC2)<(m12-(10*s12)) | as.numeric(ADRN_only$PC2)>(m12+(10*s12))] <- 1
ADRN_only$PC2_outlier[as.numeric(ADRN_only$PC2)>=(m12-(10*s12)) & as.numeric(ADRN_only$PC2)<=(m12+(10*s12))] <- 0
ADRN_only$PC3_outlier[as.numeric(ADRN_only$PC3)<(m13-(6*s13)) | as.numeric(ADRN_only$PC3)>(m13+(6*s13))] <- 1
ADRN_only$PC3_outlier[as.numeric(ADRN_only$PC3)>=(m13-(6*s13)) & as.numeric(ADRN_only$PC3)<=(m13+(6*s13))] <- 0
ADRN_only$PC4_outlier[as.numeric(ADRN_only$PC4)<(m14-(6*s14)) | as.numeric(ADRN_only$PC4)>(m14+(6*s14))] <- 1
ADRN_only$PC4_outlier[as.numeric(ADRN_only$PC4)>=(m14-(6*s14)) & as.numeric(ADRN_only$PC4)<=(m14+(6*s14))] <- 0
ADRN_only$PC5_outlier[as.numeric(ADRN_only$PC5)<(m15-(6*s15)) | as.numeric(ADRN_only$PC5)>(m15+(6*s15))] <- 1
ADRN_only$PC5_outlier[as.numeric(ADRN_only$PC5)>=(m15-(6*s15)) & as.numeric(ADRN_only$PC5)<=(m15+(6*s15))] <- 0

ADRN_outliers = subset(ADRN_only,PC1_outlier==1 | PC2_outlier==1 | PC3_outlier==1 | PC4_outlier==1 | PC5_outlier==1)
write.table(ADRN_outliers,"ADRN_outliers.txt",row.names=F,col.names=T,quote=F,sep="\t")

ADRN_outliers <- read.table("ADRN_outliers.txt",header=T,all=T,sep="\t")
ADRN_outliers2 = subset(ADRN_outliers,select=c("PATIENT","AFFSTAT","PC1_outlier","PC2_outlier"))

ADRN_outliers2$EXTRACOL <- 1
ADRN2 <- merge(ADRN_only,ADRN_outliers2,all=T)
ADRNfin <- ADRN2[is.na(ADRN2$EXTRACOL),]
ADRNfin$EXTRACOL <- NULL
ADRN_final = subset(ADRNfin,select=c("PATIENT","PC1","PC2","PC3","PC4","PC5"))
write.table(ADRN_final,"ADRN_WGS_779.txt",row.names=F,col.names=T,quote=F,sep="\t")

#Create table of outliers
library(grid)
grid.newpage()
library(gridExtra)
ADRN_theme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = 0.75)),
                                        colhead = list(fg_params=list(cex = 0.75)),
                                        rowhead = list(fg_params=list(cex = 0.75)))
grid.table(ADRN_outliers2,rows=NULL,theme=ADRN_theme)
dev.off()
