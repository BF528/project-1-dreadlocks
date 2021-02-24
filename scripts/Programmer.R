# loading packages 
library("affy")
library("affyPLM")
library("sva")
library("AnnotationDbi")
library("hgu133plus2.db")
library("ggplot2")
library("cowplot")
path_out = '/projectnb/bf528/users/dreadlocks/project_1/results/'
load("/projectnb/bf528/users/dreadlocks/project_1/scripts/sarahs_environment.RData")

# normalization
data <- ReadAffy(celfile.path = "/projectnb/bf528/users/dreadlocks/project_1/samples") # read CEL files
norm_eset <- rma(data) # create an ExpressionSet
fdata <- fitPLM(data, normalize=TRUE, background=TRUE) # This function converts an AffyBatch into an PLMset by fitting a specified robust linear model to the probe level data

# visualization
RLE(fdata, main="RLE for dataset", xaxt='n')
NUSE(fdata,main="NUSE for dataset", xaxt='n')
statsRLE <- RLE(fdata, type="stats")
statsNUSE <- NUSE(fdata, type="stats")
medianRLE <- statsRLE["median",]
medianNUSE <- statsNUSE["median",]
par(mfrow=c(1,2))
hist(medianRLE, xlab="Median RLE", main="Histogram of Median RLE")
hist(medianNUSE, xlab="Median NUSE", main="Histogram of Median NUSE")

# correct for batch effects
metadata <- read.csv("/projectnb/bf528/users/dreadlocks/project_1/proj_metadata.csv", as.is=T)
mod <- model.matrix(~as.factor(normalizationcombatmod), data = metadata)
# batch corrected expression
batch_corrected <-ComBat(exprs(norm_eset), 
                         batch=metadata$normalizationcombatbatch, 
                         mod=mod)
# remove outliers
drop <- c(names(medianNUSE[which(medianNUSE > 1.05)]), names(medianRLE[which(medianRLE > 0.10)]))
filtered_batch_corrected <- batch_corrected[,!colnames(batch_corrected) %in% drop]
#Write out the expression data to a CSV file 
write.csv(filtered_batch_corrected,paste(path_out, 'expression_matrix.csv',sep = ''))

# check that the order of metadata and probe data are the same
metadata$samplename <- sapply(strsplit(meta_data$CEL_path, "/"), function(x) x[8])
#all.equal(meta_data_name, rownames(pData(data)))
filtered_metadata <- metadata[medianNUSE < 1.05 & medianRLE < 0.10,]

# scaling the transposed matrix 
sdata <- t(scale(t(filtered_batch_corrected)))
non_batch_sdata <- t(scale(t(exprs(norm_eset))))
# PCA
# need to transpose again because prcomp() takes columns as features
pca_batch_res <- prcomp(t(sdata), scale.=FALSE, center=FALSE)
pca_non_batch_res <- prcomp(t(non_batch_sdata), scale.=FALSE, center=FALSE)
# create a data frame that includes subtype info for plotting 
pca_batch_data <- data.frame(PC1=pca_batch_res$x[,"PC1"], 
                             PC2=pca_batch_res$x[,"PC2"],
                             Subtype = filtered_metadata$SixSubtypesClassification)
pca_non_batch_data <- data.frame(PC1=pca_non_batch_res$x[,"PC1"], 
                             PC2=pca_non_batch_res$x[,"PC2"],
                             Subtype = metadata$SixSubtypesClassification)

# get % variation explained by PC1 and PC2
summary(pca_batch_res)
summary(pca_non_batch_res)

# plot PCA batch corrected and not batch corrected 
pca_batch <-ggplot(data = pca_batch_data) +
  geom_point(aes(PC1, PC2, color = Subtype)) +
  xlab("PC1(11.4%)") + 
  ylab("PC2(8.4%)") +
  ggtitle("PCA With Batch Correction")
pca_non_batch <- ggplot(data = pca_non_batch_data) +
  geom_point(aes(PC1, PC2, color = Subtype)) +
  xlab("PC1(14.5%)") + 
  ylab("PC2(9.5%)") +
  ggtitle("PCA Without Batch Correction")
plot_grid(pca_batch,pca_non_batch, labels = "AUTO")


