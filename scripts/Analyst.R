# set working directory
setwd("/projectnb2/bf528/users/dreadlocks/project_1/step_4&5/")

# read the RMA normalized, ComBat adjusted gene expression values
NE = read.csv(file="/projectnb2/bf528/users/dreadlocks/project_1/results/expression_matrix.csv")

# each row represents a gene, each column represents a patient

##########------------STEP4---------------
# 4.1
# create a vector has a length of genes number
v4_1 = vector(length = length(rownames(NE)))
for(i in 1:length(rownames(NE))){
  if(sum(NE[i,2:length(colnames(NE))] > log(15, base = 2)) > (length(colnames(NE))-1)*0.2) {v4_1[i] = 1}
  else {v4_1[i] = 0}
}
# for here, 1 in v4_1 means that gene i satisfy the 4.1

# 4.2
# first, compute the median of variance of all genes
probe_variance = vector(length = length(rownames(NE)))
for(i in 1:length(rownames(NE))){
  probe_variance[i] = var(unlist(NE[i,2:length(colnames(NE))]))
}
med_var = median(probe_variance)

# compute the test statistic for each gene
test_statistic = vector(length = length(rownames(NE)))
for(i in 1:length(rownames(NE))){
  test_statistic[i] = (length(colnames(NE))-2)*(sd(unlist(NE[i,2:length(colnames(NE))])) /med_var^0.5)^2
}
# figure out which gene satisfy 4.2
# we use upper one-tailed test
v4_2 = vector(length = length(rownames(NE)))
for(i in 1:length(rownames(NE))){
  if(test_statistic[i] > qchisq(0.99, length(colnames(NE))-2)) {v4_2[i] = 1}
  else {v4_2[i] = 0}
}

#4.3
v4_3 = vector(length = length(rownames(NE)))
for(i in 1:length(rownames(NE))){
  if(sd(unlist(NE[i,2:length(colnames(NE))])) / mean((unlist(NE[i,2:length(colnames(NE))]))) > 0.186) {v4_3[i] = 1}
  else {v4_3[i] = 0}
}

#4.4
GEM4_4 = NE[which(v4_1==1 & v4_2==1 & v4_3==1),]
write.csv(GEM4_4, "/projectnb2/bf528/users/dreadlocks/project_1/results/GEM4_4.csv", row.names = FALSE)

#4.5
GEM4_2 = NE[which(v4_2==1),]
write.csv(GEM4_2, "/projectnb2/bf528/users/dreadlocks/project_1/results/GEM4_2.csv", row.names = FALSE)

rm(med_var, test_statistic, probe_variance)

##########------------STEP5---------------
#5.1
# compute the distance and plot the hc
distmat <- dist(t(GEM4_4[,2:length(colnames(NE))]), method = "euclidean")
hclust_patient <- hclust(distmat, method = "average")
plot(hclust_patient, cex = 0.6, hang = -1, main = "Hierarchical clustering of patients")

#5.2
# cut the hc and table to see the number in each group
cut_hclust <- cutree(hclust_patient, k=2)
table(cut_hclust)

#5.3
meta_data <- read.csv("/project/bf528/project_1/doc/proj_metadata.csv")
subtypevector = vector(length = length(colnames(GEM4_4))-1)
t <- 1
for(i in 1:length(meta_data$geo_accession)){
  if(sum(grepl(meta_data$geo_accession[i], colnames(GEM4_4))) == 1) {
    subtypevector[t] <- meta_data$cit.coloncancermolecularsubtype[i]
    t <- t+1
  }
}
colorvector = ifelse(subtypevector=="C3", "red", "blue")
num_mat <- matrix(as.numeric(unlist(GEM4_4[,2:length(colnames(GEM4_4))])),nrow=nrow(GEM4_4[,2:length(colnames(GEM4_4))]))
heatmap(scale = "row", num_mat, ColSideColors=colorvector, main = 'Heatmap of the gene-expression')

#5.4
p_value <- vector(length = length(rownames(GEM4_4)))
statistic <- vector(length = length(rownames(GEM4_4)))
for(i in 1:length(rownames(GEM4_4))){
  p_value[i] <- t.test(GEM4_4[i,2:length(colnames(GEM4_4))][which(as.numeric(cut_hclust)==1)], GEM4_4[i,2:length(colnames(GEM4_4))][which(as.numeric(cut_hclust)==2)])$p.value
  statistic[i] <- t.test(GEM4_4[i,2:length(colnames(GEM4_4))][which(as.numeric(cut_hclust)==1)], GEM4_4[i,2:length(colnames(GEM4_4))][which(as.numeric(cut_hclust)==2)])$statistic
}
p_adjust <- p.adjust(p_value, method="fdr")

#create a dataframe
GEM5_4 <- data.frame(GEM4_4$X, statistic, p_value, p_adjust)
names(GEM5_4) <- c("ID", "statistic", "p_value", "p_adjust")
GEM5_4 <- GEM5_4[order(GEM5_4$p_adjust),]
rm(statistic, p_value, p_adjust)
print(sum(GEM5_4$p_adjust<0.05))

#write out the dataframe
write.csv(GEM5_4, "/projectnb2/bf528/users/dreadlocks/project_1/results/GEM5_4.csv", row.names = FALSE)

#5.5
print(GEM5_4$ID[1:10])

#5.6
p_value <- vector(length = length(rownames(GEM4_2)))
statistic <- vector(length = length(rownames(GEM4_2)))
for(i in 1:length(rownames(GEM4_2))){
  p_value[i] <- t.test(GEM4_2[i,2:length(colnames(GEM4_2))][which(subtypevector=="C3")], GEM4_2[i,2:length(colnames(GEM4_2))][which(subtypevector!="C3")])$p.value
  statistic[i] <- t.test(GEM4_2[i,2:length(colnames(GEM4_2))][which(subtypevector=="C3")], GEM4_2[i,2:length(colnames(GEM4_2))][which(subtypevector!="C3")])$statistic
}
p_adjust <- p.adjust(p_value, method="fdr")

#create a dataframe
GEM5_6 <- data.frame(GEM4_2$X, statistic, p_value, p_adjust)
names(GEM5_6) <- c("ID", "statistic", "p_value", "p_adjust")
GEM5_6 <- GEM5_6[order(GEM5_6$p_adjust),]
rm(statistic, p_value, p_adjust)

#write out the dataframe
write.csv(GEM5_6, "/projectnb2/bf528/users/dreadlocks/project_1/results/GEM5_6.csv", row.names = FALSE)





