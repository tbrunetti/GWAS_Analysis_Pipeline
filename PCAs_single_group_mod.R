#source("https://bioconductor.org/biocLite.R")
#biocLite("GENESIS")
library("GENESIS")
# loaded from genesis setup script in GWAS_analysis_pipeline.py
load('/home/brunettt/GENESIS_VIF2_MAF01/individual_GENESIS/CAAPA_MEGA_pass_QC_removed_failures_include_Peru_African_American_maf_greater_thresh_hetFiltered_dups_removed_GENESIS')
get_dataframe <- pData(scanAnnot)

population_name = 'HONDURAS_NO_BELEN_no-PCA_outliers'


pdf(paste('/home/brunettt/GENESIS_VIF2_MAF01/',population_name, "_individual_PCA_plots.pdf", sep=''),onefile=T)
par(mfrow=c(1,1))
par(mar=c(5,4.5,0,2))
par(omi=c(.5,0,.5,0))
fig=c(3,9,1,4)/10

#Plot first 200 Eigenvalues from mypcair
data <- data.frame(as.matrix(mypcair$values))
data['eigenvalues'] <- data.frame(as.matrix(mypcair$values))
data['proportion_var_explained'] <-data.frame((as.matrix(mypcair$values/mypcair$sum.values)))
data['PC'] <- data.frame((as.numeric(rownames(data))))
subset_data <- subset(data,PC<201)
subset_data['new_eigensum'] <-data.frame(colSums(subset_data)['eigenvalues'])
subset_data['subset_Prop_variance'] <- data.frame(subset_data$eigenvalues/subset_data$new_eigensum)
plot(subset_data$PC, subset_data$subset_Prop_variance, col="black", pch=20, cex=1,ylab="Proportion of Variance Explained",xlab="PCs",xaxt="n")
axis(1,at=subset_data$PC,labels=subset_data$PC)

#Plot first 20 Eigenvalues from mypcair
data <- data.frame(as.matrix(mypcair$values))
data['eigenvalues'] <- data.frame(as.matrix(mypcair$values))
data['proportion_var_explained'] <-data.frame((as.matrix(mypcair$values/mypcair$sum.values)))
data['PC'] <- data.frame((as.numeric(rownames(data))))
subset_data <- subset(data,PC<21)
subset_data['new_eigensum'] <-data.frame(colSums(subset_data)['eigenvalues'])
subset_data['subset_Prop_variance'] <- data.frame(subset_data$eigenvalues/subset_data$new_eigensum)
plot(subset_data$PC, subset_data$subset_Prop_variance, col="black", pch=20, cex=1,ylab="Proportion of Variance Explained",xlab="PCs",xaxt="n")
axis(1,at=subset_data$PC,labels=subset_data$PC)

# Assign group name
pcs_with_pops <- get_dataframe
pcs_with_pops$population <- population_name

# color order for scatter plots
cols=c("blue")

# scatter PC1 vs PC2
plot(pcs_with_pops$pc1, pcs_with_pops$pc2, col=cols,
     pch=c(17), cex=0.5,xlab="Principal Component 1",ylab="Principal Component 2")
abline(v=mean(pcs_with_pops$pc1), col='darkgrey', lty=2)
abline(v=mean(pcs_with_pops$pc1)+(1*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(1*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(2*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(2*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(3*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(3*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(4*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(4*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(5*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(5*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(6*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(6*sd(pcs_with_pops$pc1)), col='darkgrey')

abline(h=mean(pcs_with_pops$pc2), col='darkgrey', lty=2)
abline(h=mean(pcs_with_pops$pc2)+(1*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc2)-(1*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc2)+(2*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc2)-(2*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc2)+(3*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc2)-(3*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc2)+(4*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc2)-(4*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc2)+(5*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc2)-(5*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc2)+(6*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc2)-(6*sd(pcs_with_pops$pc2)), col='darkgrey')
legend("bottomright", population_name ,pch=c(21),pt.bg=cols,pt.cex=1,cex=.5)

# scatter PC1 vs PC3
plot(pcs_with_pops$pc1, pcs_with_pops$pc3, col=cols,
     pch=c(17), cex=0.5,xlab="Principal Component 1",ylab="Principal Component 3")
abline(v=mean(pcs_with_pops$pc1), col='darkgrey', lty=2)
abline(v=mean(pcs_with_pops$pc1)+(1*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(1*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(2*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(2*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(3*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(3*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(4*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(4*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(5*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(5*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(6*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(6*sd(pcs_with_pops$pc1)), col='darkgrey')

abline(h=mean(pcs_with_pops$pc3), col='darkgrey', lty=2)
abline(h=mean(pcs_with_pops$pc3)+(1*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(1*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)+(2*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(2*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)+(3*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(3*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)+(4*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(4*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)+(5*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(5*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)+(6*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(6*sd(pcs_with_pops$pc3)), col='darkgrey')
legend("topright", population_name ,pch=c(21),pt.bg=cols,pt.cex=1,cex=.5)

# scatter PC1 vs PC4
plot(pcs_with_pops$pc1, pcs_with_pops$pc4, col=cols,
     pch=c(17), cex=0.5,xlab="Principal Component 1",ylab="Principal Component 4")
abline(v=mean(pcs_with_pops$pc1), col='darkgrey', lty=2)
abline(v=mean(pcs_with_pops$pc1)+(1*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(1*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(2*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(2*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(3*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(3*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(4*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(4*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(5*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(5*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(6*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(6*sd(pcs_with_pops$pc1)), col='darkgrey')

abline(h=mean(pcs_with_pops$pc4), col='darkgrey', lty=2)
abline(h=mean(pcs_with_pops$pc4)+(1*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(1*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(2*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(2*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(3*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(3*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(4*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(4*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(5*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(5*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(6*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(6*sd(pcs_with_pops$pc4)), col='darkgrey')

legend("bottomright", population_name,pch=c(21),pt.bg=cols,pt.cex=1,cex=.5)


# scatter PC1 vs PC5
plot(pcs_with_pops$pc1, pcs_with_pops$pc5, col=cols,
     pch=c(17), cex=0.5,xlab="Principal Component 1",ylab="Principal Component 5")
abline(v=mean(pcs_with_pops$pc1), col='darkgrey', lty=2)
abline(v=mean(pcs_with_pops$pc1)+(1*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(1*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(2*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(2*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(3*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(3*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(4*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(4*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(5*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(5*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)+(6*sd(pcs_with_pops$pc1)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc1)-(6*sd(pcs_with_pops$pc1)), col='darkgrey')

abline(h=mean(pcs_with_pops$pc5), col='darkgrey', lty=2)
abline(h=mean(pcs_with_pops$pc5)+(1*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(1*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(2*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(2*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(3*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(3*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(4*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(4*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(5*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(5*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(6*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(6*sd(pcs_with_pops$pc5)), col='darkgrey')
legend("bottomright", population_name,pch=c(21),pt.bg=cols,pt.cex=1,cex=.5)

# scatter PC2 vs PC3
plot(pcs_with_pops$pc2, pcs_with_pops$pc3, col=cols,
     pch=c(17), cex=0.5,xlab="Principal Component 2",ylab="Principal Component 3")
abline(v=mean(pcs_with_pops$pc2), col='darkgrey', lty=2)
abline(v=mean(pcs_with_pops$pc2)+(1*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(1*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(2*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(2*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(3*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(3*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(4*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(4*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(5*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(5*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(6*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(6*sd(pcs_with_pops$pc2)), col='darkgrey')

abline(h=mean(pcs_with_pops$pc3), col='darkgrey', lty=2)
abline(h=mean(pcs_with_pops$pc3)+(1*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(1*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)+(2*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(2*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)+(3*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(3*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)+(4*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(4*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)+(5*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(5*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)+(6*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc3)-(6*sd(pcs_with_pops$pc3)), col='darkgrey')

legend("topleft", population_name,pch=c(21),pt.bg=cols,pt.cex=1,cex=.5)


# scatter PC2 vs PC4
plot(pcs_with_pops$pc2, pcs_with_pops$pc4, col=cols,
     pch=c(17), cex=0.5,xlab="Principal Component 2",ylab="Principal Component 4")
abline(v=mean(pcs_with_pops$pc2), col='darkgrey', lty=2)
abline(v=mean(pcs_with_pops$pc2)+(1*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(1*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(2*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(2*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(3*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(3*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(4*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(4*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(5*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(5*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(6*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(6*sd(pcs_with_pops$pc2)), col='darkgrey')

abline(h=mean(pcs_with_pops$pc4), col='darkgrey', lty=2)
abline(h=mean(pcs_with_pops$pc4)+(1*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(1*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(2*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(2*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(3*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(3*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(4*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(4*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(5*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(5*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(6*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(6*sd(pcs_with_pops$pc4)), col='darkgrey')
legend("bottomleft", population_name,pch=c(21),pt.bg=cols,pt.cex=1,cex=.5)


# scatter PC2 vs PC5
plot(pcs_with_pops$pc2, pcs_with_pops$pc5, col=cols,
     pch=c(17), cex=0.5,xlab="Principal Component 2",ylab="Principal Component 5")
abline(v=mean(pcs_with_pops$pc2), col='darkgrey', lty=2)
abline(v=mean(pcs_with_pops$pc2)+(1*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(1*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(2*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(2*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(3*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(3*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(4*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(4*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(5*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(5*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)+(6*sd(pcs_with_pops$pc2)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc2)-(6*sd(pcs_with_pops$pc2)), col='darkgrey')

abline(h=mean(pcs_with_pops$pc5), col='darkgrey', lty=2)
abline(h=mean(pcs_with_pops$pc5)+(1*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(1*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(2*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(2*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(3*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(3*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(4*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(4*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(5*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(5*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(6*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(6*sd(pcs_with_pops$pc5)), col='darkgrey')
legend("bottomleft", population_name,pch=c(21),pt.bg=cols,pt.cex=1,cex=.5)


# scatter PC3 vs PC4
plot(pcs_with_pops$pc3, pcs_with_pops$pc4, col=cols,
     pch=c(17), cex=0.5,xlab="Principal Component 3",ylab="Principal Component 4")
abline(v=mean(pcs_with_pops$pc3), col='darkgrey', lty=2)
abline(v=mean(pcs_with_pops$pc3)+(1*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(1*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)+(2*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(2*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)+(3*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(3*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)+(4*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(4*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)+(5*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(5*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)+(6*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(6*sd(pcs_with_pops$pc3)), col='darkgrey')

abline(h=mean(pcs_with_pops$pc4), col='darkgrey', lty=2)
abline(h=mean(pcs_with_pops$pc4)+(1*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(1*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(2*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(2*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(3*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(3*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(4*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(4*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(5*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(5*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)+(6*sd(pcs_with_pops$pc4)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc4)-(6*sd(pcs_with_pops$pc4)), col='darkgrey')
legend("bottomleft", population_name, pch=c(21),pt.bg=cols,pt.cex=1,cex=.5)


# scatter PC3 vs PC5
plot(pcs_with_pops$pc3, pcs_with_pops$pc5, col=cols,
     pch=c(17), cex=0.5,xlab="Principal Component 3",ylab="Principal Component 5")
abline(v=mean(pcs_with_pops$pc3), col='darkgrey', lty=2)
abline(v=mean(pcs_with_pops$pc3)+(1*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(1*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)+(2*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(2*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)+(3*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(3*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)+(4*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(4*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)+(5*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(5*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)+(6*sd(pcs_with_pops$pc3)), col='darkgrey')
abline(v=mean(pcs_with_pops$pc3)-(6*sd(pcs_with_pops$pc3)), col='darkgrey')

abline(h=mean(pcs_with_pops$pc5), col='darkgrey', lty=2)
abline(h=mean(pcs_with_pops$pc5)+(1*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(1*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(2*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(2*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(3*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(3*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(4*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(4*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(5*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(5*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)+(6*sd(pcs_with_pops$pc5)), col='darkgrey')
abline(h=mean(pcs_with_pops$pc5)-(6*sd(pcs_with_pops$pc5)), col='darkgrey')
legend("bottomleft", population_name,pch=c(21),pt.bg=cols,pt.cex=1,cex=.5)




#Boxplots centered on mean African ancestry
PCA_pops_AFR = subset(pcs_with_pops,population==population_name)


#boxplot for PC1
m11=mean(PCA_pops_AFR$pc1)
s11=sd(PCA_pops_AFR$pc1)
boxplot(pc1 ~ population, data = pcs_with_pops,ylab = "PC1", xlab=paste("PC1, +/- 6 SD on",population_name, sep=' '),cex=0.5,outline=F)
stripchart(pc1 ~ population, data = pcs_with_pops,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m11, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m11-(6*s11), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m11+(6*s11), col = "black", lty=2, lwd=1,xpd=F)

#boxplot for PC2
m12=mean(PCA_pops_AFR$pc2)
s12=sd(PCA_pops_AFR$pc2)
boxplot(pc2 ~ population, data = pcs_with_pops,ylab = "PC2", xlab=paste("PC2, +/- 6 SD on", population_name, sep=' '),cex=0.5,outline=F)
stripchart(pc2 ~ population, data = pcs_with_pops,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m12, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m12-(6*s12), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m12+(6*s12), col = "black", lty=2, lwd=1,xpd=F)

#boxplot for PC3
m13=mean(PCA_pops_AFR$pc3)
s13=sd(PCA_pops_AFR$pc3)
boxplot(pc3 ~ population, data = pcs_with_pops,ylab = "PC3", xlab=paste("PC3, +/- 6 SD on", population_name, sep=' '),cex=0.5,outline=F)
stripchart(pc3 ~ population, data = pcs_with_pops,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m13, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m13-(6*s13), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m13+(6*s13), col = "black", lty=2, lwd=1,xpd=F)

#boxplot for PC4
m14=mean(PCA_pops_AFR$pc4)
s14=sd(PCA_pops_AFR$pc4)
boxplot(pc4 ~ population, data = pcs_with_pops,ylab = "PC4", xlab=paste("PC4, +/- 6 SD on", population_name, sep=' '),cex=0.5,outline=F)
stripchart(pc4 ~ population, data = pcs_with_pops,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m14, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m14-(6*s14), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m14+(6*s14), col = "black", lty=2, lwd=1,xpd=F)

#boxplot for PC5
m15=mean(PCA_pops_AFR$pc5)
s15=sd(PCA_pops_AFR$pc5)
boxplot(pc5 ~ population, data = pcs_with_pops,ylab = "PC5", xlab=paste("PC5, +/- 6 SD on", population_name, sep=' '),cex=0.5,outline=F)
stripchart(pc5 ~ population, data = pcs_with_pops,vertical=T,method="jitter",add=T,pch=20,col=cols,cex=0.5)
abline(h =m15, col = "black", lty=1, lwd=2,xpd=F)
abline(h =m15-(6*s15), col = "black", lty=2, lwd=1,xpd=F)
abline(h =m15+(6*s15), col = "black", lty=2, lwd=1,xpd=F)


#Determine outliers
BASS_only = subset(pcs_with_pops,population==population_name)
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
write.table(BASS_outliers,paste(population_name,"_outliers_based_on_individual_pcs.txt", sep=''),row.names=F,col.names=T,quote=F,sep="\t")

BASS_outliers <- read.table(paste(population_name,"_outliers_based_on_individual_pcs.txt", sep=''),header=T,all=T,sep="\t")
BASS_outliers2 = subset(BASS_outliers,select=c("scanID","pheno","PC1_outlier","PC2_outlier"))

BASS_outliers2$EXTRACOL <- 1
BASS2 <- merge(BASS_only,BASS_outliers2,all=T)
BASSfin <- BASS2[is.na(BASS2$EXTRACOL),]
BASSfin$EXTRACOL <- NULL
BASS_final = subset(BASSfin,select=c("scanID","pc1","pc2","pc3","pc4","pc5"))
write.table(BASS_final,paste(population_name,"_PCs_No_Outliers_based_on_individual_pcs.txt", sep=''),row.names=F,col.names=T,quote=F,sep="\t")


#Create table of outliers
#library(grid)
#grid.newpage()
#library(gridExtra)
#ADRN_theme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = 0.2)),
#                                        colhead = list(fg_params=list(cex = 0.2)),
#                                        rowhead = list(fg_params=list(cex = 0.2)))
#grid.table(BASS_outliers2,rows=NULL,theme=ADRN_theme)
dev.off()

#------------------------------------ONLY FOR FINDING PCs ASSOCIATED WITH PHENO-----------------------------------------
# this is updated to a logistic regression model and not a linear model using GENESIS
# loaded from genesis setup script in GWAS_analysis_pipeline.py; loading automatically creates scanAnnot variable object
load('/home/brunettt/GENESIS_VIF2_MAF01/GENESIS_changed_name_to_match_vcf/CAAPA_MEGA_pass_QC_removed_failures_include_Peru_African_American_maf_greater_thresh_hetFiltered_dups_removed_GENESIS') 
population_name = 'JAMAICA'
all_pvalues = list()
#df=as.data.frame(mypcair$vectors)
#merged_df <- merge(df, df_original, by=0, all = TRUE)
#no_sample_names = subset(merged_df, select = -c(Row.names))
#no_sample_names$ASTHMA=no_sample_names$pheno
pcs_associated_with_asthma=list()
for (i in 1:20) {
  pc_name <- paste0("pc",i)
  nullmod <- fitNullMM(scanData=scanAnnot, 
                       outcome='pheno',
                       covars=c(pc_name), 
                       covMatList=covMatList,
                       family=binomial(link='logit'))
  pvalue <- nullmod$fixef[2,"pval"]
  all_pvalues <- c(all_pvalues, pvalue)
  if (pvalue<0.05) {
    pcs_associated_with_asthma <- c(pcs_associated_with_asthma, pc_name)
  }
}
  
#df2=do.call("rbind",l1)
#write.table(df2,file=paste('/home/brunettt/GENESIS_VIF2_MAF01/logistic_regressions/',population_name,"_LM_pvalues_PC_ASTHMA.txt", sep=''),sep="\t",col.names=F,row.names=F,quote=F)



