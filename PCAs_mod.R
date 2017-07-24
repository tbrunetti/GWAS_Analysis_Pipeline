load('/home/tonya/Desktop/CAAPA_MEGA_passed_QC_final_converted_Honduran_maf_greater_thresh_hetFiltered_dups_removed_thousGen_GENESIS')
get_dataframe <- pData(scanAnnot)
write.table(get_dataframe, file='/home/tonya/Desktop/get_dataframe.txt', sep = '\t')


pcs_with_pops <- read.table('/home/tonya/Desktop/merged_pop_Honduran_CAAPA.txt')

pdf("HONDURAN_PCA_plots.pdf",onefile=T)
par(mfrow=c(1,1))
par(mar=c(5,4.5,0,2))
par(omi=c(.5,0,.5,0))
fig=c(3,9,1,4)/10

#Plot first N Eigenvalues from mypcair
data <- data.frame(as.matrix(mypcair$values))
data['eigenvalues'] <- data.frame(as.matrix(mypcair$values))
data['proportion_var_explained'] <-data.frame((as.matrix(mypcair$values/mypcair$sum.values)))
data['PC'] <- data.frame((as.numeric(rownames(data))))
subset_data <- subset(data,PC<201)
subset_data['new_eigensum'] <-data.frame(colSums(subset_data)['eigenvalues'])
subset_data['subset_Prop_variance'] <- data.frame(subset_data$eigenvalues/subset_data$new_eigensum)
plot(subset_data$PC, subset_data$subset_Prop_variance, col="black", pch=20, cex=1,ylab="Proportion of Variance Explained",xlab="",xaxt="n")
axis(1,at=subset_data$PC,labels=subset_data$PC)

# Assign experimental vs TGP
pcs_with_pops$group[pcs_with_pops$population=="HONDURAN"] <- "HONDURAN"
pcs_with_pops$group[pcs_with_pops$population!="HONDURAN"] <- "TGP"
pcs_with_pops$group <- as.factor(pcs_with_pops$group)

# scatter PC1 vs PC2
cols=c("blue", "darkgreen", "orange", "magenta", "pink", "red")
plot(pcs_with_pops$pc1, pcs_with_pops$pc2, col=cols[as.numeric(pcs_with_pops$population)],
     pch=c(17,20)[as.numeric(pcs_with_pops$group)], cex=0.5,xlab="Principal Component 1",ylab="Principal Component 2")
legend("topright", c("AFR","AMR","EAS","EUR","HONDURAN","SAS"),pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=1,cex=.5)
legend("bottomleft", c("HONDURAN","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)

# scatter PC1 vs PC3
cols=c("blue", "darkgreen", "orange", "magenta", "pink", "red")
plot(pcs_with_pops$pc1, pcs_with_pops$pc3, col=cols[as.numeric(pcs_with_pops$population)],
     pch=c(17,20)[as.numeric(pcs_with_pops$group)], cex=0.5,xlab="Principal Component 1",ylab="Principal Component 3")
legend("topright", c("AFR","AMR","EAS","EUR","HONDURAN","SAS"),pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=1,cex=0.5)
legend("bottomleft", c("HONDURAN","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)

# scatter PC1 vs PC4
cols=c("blue", "darkgreen", "orange", "magenta", "pink", "red")
plot(pcs_with_pops$pc1, pcs_with_pops$pc4, col=cols[as.numeric(pcs_with_pops$population)],
     pch=c(17,20)[as.numeric(pcs_with_pops$group)], cex=0.5,xlab="Principal Component 1",ylab="Principal Component 4")
legend("topright", c("AFR","AMR","EAS","EUR","HONDURAN","SAS"),pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=1,cex=0.5)
legend("bottomleft", c("HONDURAN","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)


# scatter PC1 vs PC5
cols=c("blue", "darkgreen", "orange", "magenta", "pink", "red")
plot(pcs_with_pops$pc1, pcs_with_pops$pc5, col=cols[as.numeric(pcs_with_pops$population)],
     pch=c(17,20)[as.numeric(pcs_with_pops$group)], cex=0.5,xlab="Principal Component 1",ylab="Principal Component 5")
legend("topright", c("AFR","AMR","EAS","EUR","HONDURAN","SAS"),pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=1,cex=0.5)
legend("bottomleft", c("HONDURAN","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)


# scatter PC2 vs PC3
cols=c("blue", "darkgreen", "orange", "magenta", "pink", "red")
plot(pcs_with_pops$pc2, pcs_with_pops$pc3, col=cols[as.numeric(pcs_with_pops$population)],
     pch=c(17,20)[as.numeric(pcs_with_pops$group)], cex=0.5,xlab="Principal Component 2",ylab="Principal Component 3")
legend("topright", c("AFR","AMR","EAS","EUR","HONDURAN","SAS"),pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=1,cex=0.5)
legend("bottomleft", c("HONDURAN","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)

# scatter PC2 vs PC4
cols=c("blue", "darkgreen", "orange", "magenta", "pink", "red")
plot(pcs_with_pops$pc2, pcs_with_pops$pc4, col=cols[as.numeric(pcs_with_pops$population)],
     pch=c(17,20)[as.numeric(pcs_with_pops$group)], cex=0.5,xlab="Principal Component 2",ylab="Principal Component 4")
legend("topright", c("AFR","AMR","EAS","EUR","HONDURAN","SAS"),pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=1,cex=0.5)
legend("bottomleft", c("HONDURAN","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)

# scatter PC2 vs PC5
cols=c("blue", "darkgreen", "orange", "magenta", "pink", "red")
plot(pcs_with_pops$pc2, pcs_with_pops$pc5, col=cols[as.numeric(pcs_with_pops$population)],
     pch=c(17,20)[as.numeric(pcs_with_pops$group)], cex=0.5,xlab="Principal Component 2",ylab="Principal Component 5")
legend("topright", c("AFR","AMR","EAS","EUR","HONDURAN","SAS"),pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=1,cex=0.5)
legend("bottomleft", c("HONDURAN","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)

# scatter PC3 vs PC4
cols=c("blue", "darkgreen", "orange", "magenta", "pink", "red")
plot(pcs_with_pops$pc3, pcs_with_pops$pc4, col=cols[as.numeric(pcs_with_pops$population)],
     pch=c(17,20)[as.numeric(pcs_with_pops$group)], cex=0.5,xlab="Principal Component 3",ylab="Principal Component 4")
legend("topright", c("AFR","AMR","EAS","EUR","HONDURAN","SAS"),pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=1,cex=0.5)
legend("bottomleft", c("HONDURAN","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)


# scatter PC3 vs PC5
cols=c("blue", "darkgreen", "orange", "magenta", "pink", "red")
plot(pcs_with_pops$pc3, pcs_with_pops$pc5, col=cols[as.numeric(pcs_with_pops$population)],
     pch=c(17,20)[as.numeric(pcs_with_pops$group)], cex=0.5,xlab="Principal Component 3",ylab="Principal Component 5")
legend("topright", c("AFR","AMR","EAS","EUR","HONDURAN","SAS"),pch=c(21,21,21,21,21,21),pt.bg=cols,pt.cex=1,cex=0.5)
legend("bottomleft", c("HONDURAN","TGP"),pch=c(17,20),pt.cex=0.5,cex=0.5)

#Boxplots centered on mean African ancestry
PCA_pops_AFR = subset(pcs_with_pops,population=="AFR")
#boxplot for PC1
m11=mean(PCA_pops_AFR$pc1)
s11=sd(PCA_pops_AFR$pc1)
boxplot(pc1 ~ population, data = pcs_with_pops,ylab = "PC1", xlab="PC1, +/- 6 SD on AFR",cex=0.5,outline=F)
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
BASS_only = subset(pcs_with_pops,population=="HONDURAN")
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
write.table(BASS_outliers,"HONDURAN_outliers.txt",row.names=F,col.names=T,quote=F,sep="\t")

BASS_outliers <- read.table("HONDURAN_outliers.txt",header=T,all=T,sep="\t")
BASS_outliers2 = subset(BASS_outliers,select=c("scanID","pheno","PC1_outlier","PC2_outlier"))

BASS_outliers2$EXTRACOL <- 1
BASS2 <- merge(BASS_only,BASS_outliers2,all=T)
BASSfin <- BASS2[is.na(BASS2$EXTRACOL),]
BASSfin$EXTRACOL <- NULL
BASS_final = subset(BASSfin,select=c("scanID","pc1","pc2","pc3","pc4","pc5"))
write.table(BASS_final,"HONDURAN_PCs_No_Outliers.txt",row.names=F,col.names=T,quote=F,sep="\t")


#Create table of outliers
#library(grid)
#grid.newpage()
#library(gridExtra)
#ADRN_theme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = 0.2)),
#                                        colhead = list(fg_params=list(cex = 0.2)),
#                                        rowhead = list(fg_params=list(cex = 0.2)))
#grid.table(BASS_outliers2,rows=NULL,theme=ADRN_theme)
dev.off()
