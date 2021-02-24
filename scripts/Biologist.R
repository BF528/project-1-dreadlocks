setwd("/projectnb/bf528/users/dreadlocks/project_1/")

if(!require(c("hgu133plus2.db", "GSEABase"))) BiocManager::install(c("hgu133plus2.db", "GSEABase"))

library(hgu133plus2.db)
library(GSEABase)
library(magrittr)

t_cutoff_2_summary <- read.csv("/projectnb2/bf528/users/dreadlocks/project_1/results/GEM5_6.csv",
                               row.names = 1)

gene_symbol <- select(hgu133plus2.db, keys = rownames(t_cutoff_2_summary),
                      columns = c("PROBEID","SYMBOL"), keytype = "PROBEID")
t_cutoff_2_summary$symbol <- gene_symbol[!duplicated(gene_symbol$PROBEID), "SYMBOL"]
all.equal(rownames(t_cutoff_2_summary), gene_symbol[!duplicated(gene_symbol$PROBEID), "PROBEID"])
head(t_cutoff_2_summary)

deg_table <- rbind(t_cutoff_2_summary[order(t_cutoff_2_summary$statistic, decreasing = T),][1:10,],
                   t_cutoff_2_summary[order(t_cutoff_2_summary$statistic, decreasing = F),][1:10,])
deg_table

kegg <- getGmt("/projectnb/bf528/users/dreadlocks/project_1/results/c2.cp.kegg.v7.2.symbols.gmt")
kegg
go <- getGmt("/projectnb/bf528/users/dreadlocks/project_1/results/c5.go.v7.2.symbols.gmt")
go
hallmark <- getGmt("/projectnb/bf528/users/dreadlocks/project_1/results/h.all.v7.2.symbols.gmt")
hallmark

## gsea analysis
find_overlap <- function(geneset, geneset_collection){
  
  index <- order(abs(geneset$statistic), decreasing = T)[1:1000]
  deg <- geneset[index, "symbol"]
  non_deg <- geneset[-index, "symbol"]
  
  ingene_deg <- sum(deg %in% geneset_collection)
  ingene_non_deg <- sum(non_deg %in% geneset_collection)
  
  not_ingene_deg <- sum(! deg %in% geneset_collection)
  not_ingene_non_deg <- sum(! non_deg %in% geneset_collection)
  
  res <- fisher.test(matrix(c(ingene_deg,not_ingene_deg,ingene_non_deg,not_ingene_non_deg),nrow=2))
  return(c(res$estimate,res$p.value))

}

wrapper <- function(query, genesetcollection){
  
  res <- sapply(seq_len(length(genesetcollection@.Data)), function(x){
    
    p <- find_overlap(query, geneIds(genesetcollection@.Data[[x]]))
    tmp_res <- c(setName(genesetcollection@.Data[[x]]), p)
    names(tmp_res) <- c("Term", "Odds_ratio", "p-value")
    return(tmp_res)
    
  })
  
  res <- t(res) %>% data.frame()
  res$p_adj <- p.adjust(res$p.value, method = "BH")
  res <- res[order(res$p_adj, decreasing = F),]
  
  return(res)
  
}

go_up_res <- wrapper(t_cutoff_2_summary[t_cutoff_2_summary$statistic > 0,],
                     go)
go_down_res <- wrapper(t_cutoff_2_summary[t_cutoff_2_summary$statistic < 0,],
                       go)

kegg_up_res <- wrapper(t_cutoff_2_summary[t_cutoff_2_summary$statistic > 0,],
                       kegg)
kegg_down_res <- wrapper(t_cutoff_2_summary[t_cutoff_2_summary$statistic < 0,],
                         kegg)

hallmark_up_res <- wrapper(t_cutoff_2_summary[t_cutoff_2_summary$statistic > 0,],
                           hallmark)
hallmark_down_res <- wrapper(t_cutoff_2_summary[t_cutoff_2_summary$statistic < 0,],
                             hallmark)
## import
write.csv(go_up_res, file = "/projectnb/bf528/users/dreadlocks/project_1/results/go_up_res.csv")
write.csv(go_down_res, file = "/projectnb/bf528/users/dreadlocks/project_1/results/go_down_res.csv")
write.csv(kegg_up_res, file = "/projectnb/bf528/users/dreadlocks/project_1/results/kegg_up_res.csv")
write.csv(kegg_down_res, file = "/projectnb/bf528/users/dreadlocks/project_1/results/kegg_down_res.csv")
write.csv(hallmark_up_res, file = "/projectnb/bf528/users/dreadlocks/project_1/results/hallmark_up_res.csv")
write.csv(hallmark_down_res, file = "/projectnb/bf528/users/dreadlocks/project_1/results/hallmark_down_res.csv")
